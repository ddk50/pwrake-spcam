#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "imc.h"
#include "stat.h"
#include "getargs.h"

/* only for new suprime */
#define NREGION (4)

int overscansub_y(int npx, int npy,
		  float *f, float *g,
		  int efminy, int efmaxy,
		  int osminy, int osmaxy,
		  int xmin, int xmax)
{
  int i,j,n;
  float *dat;
  float v;

  dat=(float*)malloc((osmaxy-osminy+1)*sizeof(float));

  for(i=xmin;i<=xmax;i++)
    {
      n=0;
      for(j=osminy;j<=osmaxy;j++)
	{
	  /* in fact, pixignr check is required, but not impl. */
	  dat[n]=f[i+j*npx];
	  n++;
	}
      v=floatmedian(n,dat);

      for(j=efminy;j<=efmaxy;j++)
	g[i+j*npx]=f[i+j*npx]-v;
      for(j=osminy;j<=osmaxy;j++)
	g[i+j*npx]=f[i+j*npx]-v;
    }
  free(dat);
  return 0;
}


int overscansub_x(int npx, int npy,
		  float *f, float *g,
		  int efminx, int efmaxx,
		  int osminx, int osmaxx,
		  int ymin, int ymax)
{
  int i,j,n;
  float *dat;
  float v;

  dat=(float*)malloc((osmaxx-osminx+1)*sizeof(float));

  for(j=ymin;j<=ymax;j++)
    {
      n=0;
      for(i=osminx;i<=osmaxx;i++)
	{
	  /* in fact, pixignr check is required, but not impl.*/
	  dat[n]=f[i+j*npx];
	  n++;
	}
      v=floatmedian(n,dat);
      
      for(i=efminx;i<=efmaxx;i++)
	g[i+j*npx]=f[i+j*npx]-v;
      for(i=osminx;i<=osmaxx;i++)
	g[i+j*npx]=f[i+j*npx]-v;
    }
  free(dat);
  return 0;
}

int main(int argc,char *argv[])
{
  struct imh imhin={""},imhout={""};
  struct icom icomin={0},icomout={0};

  char fnamin[100]="";
  char fnamout[100]="";
  char tmp[100]="";

  FILE *fp,*fp2;
  float *f,*g,*h;
  int i;
  int xmin=-1,xmax=-1,ymin=-1,ymax=-1;
  int ymin2,ymax2;
  int offst;
  float bzero=0,bscale=1.0,bzero_new=FLT_MIN,bscale_new=FLT_MIN;
  float blank_new=FLT_MIN;
  float pixignr=INT_MIN;
  int fitsdtype=-1;

  int npx,npy,npx2,npy2,npxos,npyos,npxef[NREGION],npyef[NREGION];
  int npxef0=0;
  int xflip,yflip;
  int region;
  char regid[2],key[10],comment[80];

  int efminx[NREGION],efmaxx[NREGION],efminy[NREGION],efmaxy[NREGION];
  int osminx[NREGION],osmaxx[NREGION],osminy[NREGION],osmaxy[NREGION];

  int j=0;
  char dtype[BUFSIZ]="";

  int notrim=0; 
  int no_y=0;

  /* test */
  int flag_norm=0;

  float gain[NREGION]={0};


  int fit_param=-1;


  getargopt opts[20];
  char *files[3]={NULL};
  int nopt=0;
  int helpflag;

  files[0]=fnamin;
  files[1]=fnamout;


  /* only trim version is supported now */
  setopts(&opts[nopt++],"-notrim", OPTTYP_FLAG , &notrim,
	  "not trim off(default:trim)",1);
  setopts(&opts[nopt++],"-no-y", OPTTYP_FLAG , &no_y,
	  "not subtract y-overscan(default:subtract)",1);

  setopts(&opts[nopt++],"-bzero=", OPTTYP_FLOAT , &bzero_new,
	  "bzero");
  setopts(&opts[nopt++],"-bscale=", OPTTYP_FLOAT , &bscale_new,
	  "bscale");
  setopts(&opts[nopt++],"-dtype=", OPTTYP_STRING , dtype,
	  "datatyp(FITSFLOAT,FITSSHORT...)");
  setopts(&opts[nopt++],"-pixignr=", OPTTYP_FLOAT , &blank_new,
	  "pixignr value");
  setopts(&opts[nopt++],"",0,NULL,NULL);

  helpflag=parsearg(argc,argv,opts,files,NULL);

  if(fnamout[0]=='\0')
   {
      fprintf(stderr,"Error: No input file specified!!\n");
      helpflag=1;
    }
  if(helpflag==1)
    {
      print_help("Usage: osmed4 <options> [filein] [fileout]",
		 opts,
		 "");
      exit(-1);
    }

  /*
   * ... Open Input
   */

  (void)printf("\n Input = %s\n",fnamin );

  if( (fp = imropen ( fnamin, &imhin, &icomin )) == NULL ) 
    {
      print_help("Usage: %s <options> [filein] [fileout]",
		 opts,
		 "");
      exit(-1);
    }

  npx=imhin.npx;
  npy=imhin.npy;

  /* get information here */
  npx2=0;
  npy2=0;
  npxos=0;
  npyos=0;

  imget_fits_value(&icomin,"S_XFLIP",tmp);
  printf("debug: xflip [%s]\n",tmp);
  if (tmp[0]=='T') xflip=1; else xflip=0;

  imget_fits_value(&icomin,"S_YFLIP",tmp);
  if (tmp[0]=='T') yflip=1; else yflip=0;

  xmin=npx;
  xmax=-1;
  ymin=npy;
  ymax=-1;
  ymin2=npy;
  ymax2=-1;

  for(region=0;region<NREGION;region++)
    {
      if (xflip==0)
	sprintf(regid,"%1d",region+1);
      else
	sprintf(regid,"%1d",NREGION-region);

      sprintf(key,"S_EFMN%c1",regid[0]);
      imget_fits_int_value(&icomin,key,&(efminx[region]));
      sprintf(key,"S_EFMX%c1",regid[0]);
      imget_fits_int_value(&icomin,key,&(efmaxx[region]));
      sprintf(key,"S_EFMN%c2",regid[0]);
      imget_fits_int_value(&icomin,key,&(efminy[region]));
      sprintf(key,"S_EFMX%c2",regid[0]);
      imget_fits_int_value(&icomin,key,&(efmaxy[region]));
     
      efminx[region]--;
      efmaxx[region]--;
      efminy[region]--;
      efmaxy[region]--;

      printf("debug:EF: %d %d %d %d\n",
	     efminx[region],efmaxx[region],efminy[region],efmaxy[region]);

      if(efmaxx[region]<efminx[region]) exit(-1);
      if(efmaxy[region]<efminy[region]) exit(-1);

      npxef[region]=efmaxx[region]-efminx[region]+1;
      npyef[region]=efmaxy[region]-efminy[region]+1;    
      
      if(xmin>efminx[region]) xmin=efminx[region];
      if(xmax<efmaxx[region]) xmax=efmaxx[region];
      if(ymin>efminy[region]) ymin=efminy[region];
      if(ymax<efmaxy[region]) ymax=efmaxy[region];

      if(ymin2>efminy[region]) ymin2=efminy[region];
      if(ymax2<efmaxy[region]) ymax2=efmaxy[region];

      printf("debug: npx2=%d\n",npx2);
      npx2+=npxef[region];

      printf("debug: region %d: %d x %d \n",
	     region,npxef[region],npyef[region]);

      sprintf(key,"S_OSMN%c1",regid[0]);
      imget_fits_int_value(&icomin,key,&(osminx[region]));
      sprintf(key,"S_OSMX%c1",regid[0]);
      imget_fits_int_value(&icomin,key,&(osmaxx[region]));
      sprintf(key,"S_OSMN%c2",regid[0]);
      imget_fits_int_value(&icomin,key,&(osminy[region]));
      sprintf(key,"S_OSMX%c2",regid[0]);
      imget_fits_int_value(&icomin,key,&(osmaxy[region]));

      osminx[region]--;
      osmaxx[region]--;
      osminy[region]--;
      osmaxy[region]--;

      printf("debug:OS: %d %d %d %d\n",
	     osminx[region],osmaxx[region],osminy[region],osmaxy[region]);

      if(ymin2>osminy[region]) ymin2=osminy[region];
      if(ymax2<osmaxy[region]) ymax2=osmaxy[region];

      if (npxos<osmaxx[region]-osminx[region]+1)
	npxos=osmaxx[region]-osminx[region]+1;

      if (flag_norm==1)
	{
	  sprintf(key,"S_GAIN%c",regid[0]);
	  imget_fits_float_value(&icomin,key,&(gain[region]));
	}
      if (npyos<osmaxy[region]-osminy[region]+1)
	npyos=osmaxy[region]-osminy[region]+1;
    }
  npy2=ymax-ymin+1;

  printf("debug: y %d-%d\n",ymin,ymax);
  printf("debug: effective size: %dx%d\n",npx2,npy2);

  /* output */
  

  imh_inherit(&imhin,&icomin,&imhout,&icomout,fnamout);
  if(!notrim)
    {
      imhout.npx   = npx2;
      imhout.npy   = npy2;
      imhout.ndatx = imhout.npx;
      imhout.ndaty = imhout.npy;
    }
  else
    {
      imhout.npx   = npx;
      imhout.npy   = npy;
      imhout.ndatx = imhout.npx;
      imhout.ndaty = imhout.npy;
    }

  if(dtype[0]=='\0')
    {
      imc_fits_get_dtype( &imhout, &fitsdtype, &bzero, &bscale, &offst );
      /** force FLOAT **/
      fitsdtype=FITSFLOAT;
    }
  else
    if ((fitsdtype=imc_fits_dtype_string(dtype))==-1)
      {
	printf("\nWarning: Unknown fits dtype %s. Inheriting original.\n",dtype);
	imc_fits_get_dtype( &imhout, &fitsdtype, &bzero, &bscale, &offst );
      }

  if(bzero_new!=FLT_MIN)
    {
      bzero=bzero_new;
    }
  else
    bzero=0;

  if(bscale_new!=FLT_MIN)
    {
      bscale=bscale_new;
    }
  else
    bscale=1.0;
  
  (void)imc_fits_set_dtype( &imhout, fitsdtype, bzero, bscale );

  if(blank_new!=FLT_MIN)
    imset_pixignr(&imhout,&icomout,(int)blank_new);    
  else
    {
      pixignr=(float)imget_pixignr(&imhin);
      imset_pixignr(&imhout,&icomout,(int)pixignr);
    }

  if (fit_param>0)
    {
      imaddhist(&icomout,"OSMED5: Overscan value is smoothed");
    }
  else
    {
      imaddhist(&icomout,"OSMED5: X-Overscan median is subtracted line by line.");
      if(no_y==0)
	imaddhist(&icomout,"OSMED5: Y-Overscan median is subtracted row by row.");
    }

  if (!notrim)
    {
      imaddhistf(&icomout,"OSMED5: And trimmed");
    }
  (void)printf("\n Output = %s\n",fnamout );


  if(!notrim)
    {
      imupdate_fitsf(&icomout, "EFP-MIN1",IMC_FITS_INTFORM,
		     1,
		     "Start position of effective frame in axis-1");
      
      imupdate_fitsf(&icomout, "EFP-MIN2",IMC_FITS_INTFORM,
		     1,
		     "Start position of effective frame in axis-2");
      
      imupdate_fitsf(&icomout, "EFP-RNG1",IMC_FITS_INTFORM,
		     npx2,
		     "Range of effective frame in axis-1");
      
      imupdate_fitsf(&icomout, "EFP-RNG2",IMC_FITS_INTFORM,
		     npy2,
		     "Range of effective frame in axis-2");
      imc_shift_WCS(&icomout,(float)xmin,(float)ymin);

      /* 2008/08/28 header revision more */
      /* tenuki kimeuchi */
      /* ch1->ch4 ? ch4->ch1?? */

      /* tenuki again */
      npxef0=0;
      for (i=0;i<NREGION;i++)
	{
	  if (xflip==0)	
	    j=i+1;
	  else
	    j=NREGION-i;
	 
	  sprintf(key,"S_EFMN%d1",j);
	  sprintf(comment,"MIN pixel of x-effective range for ch%d",j);
	  imupdate_fitsf(&icomout,key,IMC_FITS_INTFORM,
			 npxef0+1,comment);
	  sprintf(key,"S_EFMX%d1",j);
	  sprintf(comment,"MAX pixel of x-effective range for ch%d",j);
	  imupdate_fitsf(&icomout,key,IMC_FITS_INTFORM,
			 npxef0+npxef[i],comment);

	  npxef0+=npxef[i];

	  sprintf(key,"S_EFMN%d2",j);
	  sprintf(comment,"MIN pixel of y-effective range for ch%d",j);
	  imupdate_fitsf(&icomout,key,IMC_FITS_INTFORM,
			 1,comment);
	  sprintf(key,"S_EFMX%d2",j);
	  sprintf(comment,"MAX pixel of y-effective range for ch%d",j);
	  imupdate_fitsf(&icomout,key,IMC_FITS_INTFORM,
			 npy2,comment);
	}
    }

  f=(float*)malloc(npx*npy*sizeof(float));  
  imrall_tor(&imhin, fp, f, npx*npy );
  imclose(fp,&imhin,&icomin);

  g=(float*)calloc(npx*npy,sizeof(float));  

  /* overscan estimation and subtract */
  for(region=0;region<NREGION;region++)
    {
      overscansub_x(npx,npy,f,g,
		    efminx[region],efmaxx[region],
		    osminx[region],osmaxx[region],
		    ymin2,ymax2);
      if(no_y==0)
	{
	  /* 2008/07/26 Overscan Y added*/
	  overscansub_y(npx,npy,g,f,
			efminy[region],efmaxy[region],
			osminy[region],osmaxy[region],
			efminx[region],efmaxx[region]);
	}
    }
  if(no_y!=0)
    memcpy(f,g,npx*npy*sizeof(float));


  if( (fp2 = imwopen( fnamout, &imhout, &icomout )) == NULL )
    {
      print_help("Usage: %s <options> [filein] [fileout]",
		 opts,
		 "");
      exit(-1);
    }

  if(notrim)
    {
      imwall_rto( &imhout, fp2, h);
    }
  else
    {
      /* copy effective regions to g */
      i=0;
      for(region=0;region<NREGION;region++)
	{
	  for(j=ymin;j<=ymax;j++)
	    {
	      memcpy(g+i+(j-ymin)*npx2,f+efminx[region]+j*npx,
		     (efmaxx[region]-efminx[region]+1)*sizeof(float));
	    }
	  i+=(efmaxx[region]-efminx[region]+1);
	}
      imwall_rto( &imhout, fp2, g);	  
    }

  (void) imclose(fp2,&imhout,&icomout);
  return 0;
}
