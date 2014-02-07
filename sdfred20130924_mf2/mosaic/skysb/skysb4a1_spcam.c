#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "imc.h"
#include "getargs.h"
#include "stat.h"
#include "sort.h"
#include "postage.h"


float pxintp(float *g,int npx,int npy,float x,float y)
{
 /* pixel interpolation */
  int i1,i2,j1,j2;
  float dx1,dx2,dy1,dy2;
  float f;
 
  if (npx<1||npy<1) return 0; /* error */

  i1=(int)floor(x);
  if(i1>=npx-1)i1=npx-2;
  if(i1<0)i1=0;
  i2=i1+1;

  j1=(int)floor(y);
  if(j1>=npy-1)j1=npy-2;
  if(j1<0)j1=0;
  j2=j1+1;

  dx1=x-(float)i1;
  dy1=y-(float)j1;
  dx2=1.-dx1;
  dy2=1.-dy1;

  /* 2000/09/13 i2,j2 check needed here!! */
  if(i2>npx-1) 
    {
      i2=i1;
      dx1=1; 
      dx2=0;
    }
  if(j2>npy-1) 
    {
      j2=j1;
      dy1=1; 
      dy2=0;
    }

  f=dx2*dy2*g[i1+npx*j1]+dx1*dy2*g[i2+npx*j1]+
    dx2*dy1*g[i1+npx*j2]+dx1*dy1*g[i2+npx*j2];
  return f;
}

float sorted_nth(int ndat,float *dat,float n)
{
  /* dat is already sorted */
  float nf;
  int ni;

  if(n>=ndat-1) return dat[ndat-1];
  if (ndat==1 || n<0) return dat[0];
   
  ni=(int)floor(n);
  nf=n-(float)ni;
   
  return dat[ni]*(1.-nf)+nf*dat[ni];
}

float sorted_mad(int ndat,float *dat)
{
  /* not MAD infact, but good estimation */
  float med,q1,q3;
  float imed,i1;
  float df;

  imed=0.5*(float)(ndat-1);
  med=sorted_nth(ndat,dat,imed);
  df=imed*0.5;
  
  i1=0.25*(float)(ndat-1);
  q1=sorted_nth(ndat,dat,i1);
  q3=sorted_nth(ndat,dat,i1+imed);

  while(df>0.25)
    {
      if(q3+q1>2.*med)
	i1-=(float)df;
      else
	i1+=(float)df;
      df=ceil(df*2)*0.25;
      q1=sorted_nth(ndat,dat,i1);
      q3=sorted_nth(ndat,dat,i1+imed);
  }

  return 0.5*(q3-q1);
}

float skydet(int ndat0, float *g, float *skysigma, float pixignr)
{
  float median0,sigma0;
  float median,sigma;
  int i,ndat,imin,imax;
  float sky=0,ssigma=0;

  ndat=ndat0;
  imin=0;
  imax=ndat-1;
  /* sort */
  mos_heapsort(ndat,g);
  
  
  /* pixignr removal */
  for(i=0;i<ndat;i++)
    {
      if(g[i]!=pixignr) break;
    }
  imin=i;
  for(i=ndat-1;i>=0;i--)
    {
      if(g[i]!=pixignr) break;
    }
  imax=i;
  ndat=imax-imin+1;

  if (ndat<10) return 0;
  /* else */
  memmove(g,g+imin,ndat*sizeof(float));

  median0=sorted_nth(ndat,g,0.5*(float)(ndat-1));
  sigma0=sorted_mad(ndat,g)*MAD2SIGMA;

  median=median0;
  sigma=sigma0;

  while(1)
    {
      /* clip */
      for(i=0;i<ndat;i++)
	{
	  if(g[i]>median-3.*sigma) break;
	}
      imin=i;
      for(i=ndat-1;i>=0;i--)
	{
	  if(g[i]<median+3.*sigma) break;
	}
      imax=i;
      
      if(imin==0&&imax==ndat-1) break;
      ndat=imax-imin+1;
      if(ndat<10) break;

      memmove(g,g+imin,ndat*sizeof(float));

      median=sorted_nth(ndat,g,0.5*(float)(ndat-1));
      sigma=sorted_mad(ndat,g)*MAD2SIGMA;

    }


  if(ndat>0)
    {
      floatmeanrms(ndat,g,&sky,&ssigma);
      
      if(sigma/sigma0<0.8)
	sky=2.5*median-1.5*sky;
      
      if(skysigma!=NULL)
	{
	  /* lower MAD */
	  *skysigma=(sorted_nth(ndat,g,0.5*(float)(ndat-1))-
		     sorted_nth(ndat,g,0.25*(float)(ndat-1)))*MAD2SIGMA;
	}
    }
  else
    {
      sky=0;
      *skysigma=0;
    }

  return sky;
}

int skysub(float *g,int npx,int npy,
	    float *sky,
	    int nmesh_x, int nmesh_y,
	    int meshsiz_x,int meshsiz_y,
	    float pixignr)
{
  int x,y;
  float xx,yy,s;

  for(y=0;y<npy;y++)
    {
      if(nmesh_y==1)
	yy=0;
      else
	{
	  if (y<meshsiz_y/2) yy=0;
	  else if (y>npy-meshsiz_y/2) yy=nmesh_y-1;
	  else 
	    {
	      yy=(float)(y-(meshsiz_y/2))/(float)(npy-meshsiz_y)
		*(float)(nmesh_y-1);
	    }
	}      
      for(x=0;x<npx;x++)
	{
	  if (g[x+npx*y]!=pixignr)
	    {
	      if(nmesh_x==1)
		xx=0;
	      else
		{
		  if (x<meshsiz_x/2) xx=0;
		  else if (x>npx-meshsiz_x/2) xx=nmesh_x-1;
		  else 
		    {
		      xx=(float)(x-(meshsiz_x/2))/(float)(npx-meshsiz_x)
			*(float)(nmesh_x-1);
		    }
		} 
	      /* interpolate */
	      s=pxintp(sky,nmesh_x,nmesh_y,xx,yy);
	      g[x+npx*y]-=s;
	    }
	}
    }
  return 0;
}

int skypattern(float *g,int npx,int npy, 
	       int meshsiz_x,int meshsiz_y,float pixignr,
	       int nmesh_x,int nmesh_y,
	       float *skymesh, float *ssgmmesh)
{
  int ix,iy,y,x0,y0;
  float *g2;
  float dx=0.,dy=0.;

  if (nmesh_x>1)
    dx=(float)(npx-meshsiz_x)/(float)(nmesh_x-1);
  if (nmesh_y>1)
    dy=(float)(npy-meshsiz_y)/(float)(nmesh_y-1);

  g2=(float*)malloc(meshsiz_x*meshsiz_y*sizeof(float));
  
  for(iy=0;iy<nmesh_y;iy++)
    {
      y0=(int)floor(dy*(float)iy);
      for(ix=0;ix<nmesh_x;ix++)
	{
	  x0=(int)floor(dx*(float)ix);
	  for(y=0;y<meshsiz_y;y++)
	    {
	      memcpy(g2+meshsiz_x*y,
		     g+x0+(y0+y)*npx,
		     meshsiz_x*sizeof(float));		     
	    }
	  skymesh[ix+iy*nmesh_x]=skydet(meshsiz_x*meshsiz_y,g2,
					&ssgmmesh[ix+iy*nmesh_x],
					pixignr);
	}
    }

  free(g2);
  return 0;
}

int filter_skypattern(int nmesh_x,int nmesh_y,
		      float *sky,float *ssgm)
{
  int ix,iy,jx,jy,ndat,n;
  float med;
  float *sky2,s;
  float dat[9];
  /* threshold of ssgm */
  /* ssgm roughly follows chi^2~gamma*/
  /* version 0 */
  sky2=(float*)malloc(nmesh_x*nmesh_y*sizeof(float));

  ndat=nmesh_x*nmesh_y;
  med=floatmedian(ndat,ssgm);

  n=0;
  /* tenuki */
  for(iy=0;iy<nmesh_y;iy++)
    for(ix=0;ix<nmesh_x;ix++)
      {
	if (ssgm[ix+iy*nmesh_x]>med)
	  {
	    /* mask */
	    sky[ix+iy*nmesh_x]-=(ssgm[ix+iy*nmesh_x]-med);
	  }
      } 

  /* smooth */
  for(iy=0;iy<nmesh_y;iy++)
    {
      for(ix=0;ix<nmesh_x;ix++)
	{
	  n=0;
	  s=0;
	  for(jy=-1;jy<=1;jy++)
	    {
	      if(jy+iy<0||jy+iy>nmesh_y-1)
		continue;
	      for(jx=-1;jx<=1;jx++)
		{
		  if(jx+ix<0||jx+ix>nmesh_x-1)
		    continue;
		  dat[n]=sky[(jx+ix)+(jy+iy)*nmesh_x];
		  n++;
		}
	    }
	  if(n>0)
	    sky2[ix+iy*nmesh_x]=floatmedian(n,dat);
	  else
	    sky2[ix+iy*nmesh_x]=0;
	}
    }

  memcpy(sky,sky2,nmesh_x*nmesh_y*sizeof(float));
  free(sky2);

  return 0;
}

int skysub2(float *g, int npx2, int npy2, 
	    int meshsiz_x, int meshsiz_y,
	    float pixignr,
	    int refill,
	    float fillf, int fillsize_x, int fillsize_y,
	    float *medsky, float *skysigma)
{
  float *sky, *ssgm;
  int nmesh_x, nmesh_y;
  float msky,msigm;
  int i;

  /* sky estination */ 
   nmesh_x=(int)ceil(npx2/meshsiz_x*2);
   nmesh_y=(int)ceil(npy2/meshsiz_y*2);

  sky=(float*)malloc(nmesh_x*nmesh_y*sizeof(float));
  ssgm=(float*)malloc(nmesh_x*nmesh_y*sizeof(float));

  skypattern(g,npx2,npy2,meshsiz_x,meshsiz_y,pixignr,
	     nmesh_x,nmesh_y,sky,ssgm);

  filter_skypattern(nmesh_x,nmesh_y,sky,ssgm);

  msky=floatmedian(nmesh_x*nmesh_y,sky);
  msigm=floatmedian(nmesh_x*nmesh_y,ssgm);
  
  printf("debug: estimated %f %f\n",msky,msigm);
  /* interpolate */
  skysub(g,npx2,npy2,sky,nmesh_x,nmesh_y,
	 meshsiz_x,meshsiz_y,pixignr);

  free(sky);
  free(ssgm);

  if (refill)
    {
      /* again, with smaller mesh */
      if (fillf!=0)
	{
	  fillsize_x=meshsiz_x*fillf;
	  fillsize_y=meshsiz_y*fillf;
	}
      else 
	{
	  if (fillsize_x==0) fillsize_x=meshsiz_x*0.25;
	  if (fillsize_y==0) fillsize_y=meshsiz_y*0.25;
	}

      printf("debug: refill negative %d x %d\n",
	     fillsize_x,fillsize_y);

      nmesh_x=(int)ceil(npx2/fillsize_x*2)-1;
      nmesh_y=(int)ceil(npy2/fillsize_y*2)-1;
      
      sky=(float*)malloc(nmesh_x*nmesh_y*sizeof(float));
      ssgm=(float*)malloc(nmesh_x*nmesh_y*sizeof(float));
      
      skypattern(g,npx2,npy2,fillsize_x,fillsize_y,pixignr,
		 nmesh_x,nmesh_y,sky,ssgm);
      
      /* no filter, only supress + */
      for(i=0;i<nmesh_x*nmesh_y;i++)
	if(sky[i]>0) sky[i]=0;
      
      /* and fill  */
      skysub(g,npx2,npy2,sky,nmesh_x,nmesh_y,
	     fillsize_x,fillsize_y,pixignr);
      free(sky);
      free(ssgm);
    }
  
  *medsky=msky;
  *skysigma=msigm;

  return 0;
}



int main(int argc,char **argv)
{
  char fnamin[BUFSIZ]="";
  char fnamout[BUFSIZ]="";
  struct icom icomin={0},icomout={0};
  struct imh imhin={""},imhout={""};
  FILE *fpin,*fpout;

  int   npx,npy;
  int npx2,npy2;

  int meshsiz_x=0,meshsiz_y=0,meshsiz_xy=0;

  int pixignr=INT_MIN;
  int pixignr_org=INT_MIN;

  float *g,*h;
  int ndat;
  int iy;


  int xmin=-1,xmax=-1,ymin=-1,ymax=-1;


  int i;

  int refill=1;
  float fillf=0;
  int fillsize_x=0, fillsize_y=0;
  float msky,msigm;

  float bzero=FLT_MAX,bscale=FLT_MAX;
  int fitsdtype;
  char dtype[80]="";
  char key[10]="";

  /***** parse options ******/
  getargopt opts[30];
  char *files[3]={NULL};
  int nopt=0;

  files[0]=fnamin;
  files[1]=fnamout;

  nopt=0;
  
  setopts(&opts[nopt++],"-mesh=", OPTTYP_INT , &meshsiz_xy,
	  "mesh size (default: 100)");
  setopts(&opts[nopt++],"-meshx=", OPTTYP_INT , &meshsiz_x,
	  "mesh x-size (default: 100)");
  setopts(&opts[nopt++],"-meshy=", OPTTYP_INT , &meshsiz_y,
	  "mesh y-size (default: 100)");

  setopts(&opts[nopt++],"-nofill", OPTTYP_FLAG , &refill,
	  "correct negative portion",0);
  setopts(&opts[nopt++],"-filf=", OPTTYP_FLOAT , &fillf,
	  "negative check mesh relative size (default:0.25)");
  setopts(&opts[nopt++],"-fillx=", OPTTYP_INT , &fillsize_x,
	  "check negative size x (default:0.25 mesh_x)");
  setopts(&opts[nopt++],"-filly", OPTTYP_INT , &fillsize_y,
	  "check negative size y (default:0.25 mesh_y)");
  setopts(&opts[nopt++],"-bzero=", OPTTYP_FLOAT , &bzero,
	  "bzero");
  setopts(&opts[nopt++],"-bscale=", OPTTYP_FLOAT , &bscale,
	  "bscale");
  setopts(&opts[nopt++],"-dtype=", OPTTYP_STRING , dtype,
	  "datatyp(FITSFLOAT,FITSSHORT...)");
  setopts(&opts[nopt++],"-pixignr=", OPTTYP_INT , &pixignr,
	  "pixignr value");
  
  setopts(&opts[nopt++],"",0,NULL,NULL);

  if(parsearg(argc,argv,opts,files,NULL)||fnamout[0]=='\0')
    {
      print_help("Usage: skysb4 [options] filein fileout",
		 opts,"");
      exit(-1);
    }

  /**** read image ****/
  if ((fpin= imropen(fnamin,&imhin,&icomin))==NULL) 
    {
      fprintf(stderr,"File %s not found !!\n",fnamin);
      exit(-1);
    }
  
  pixignr_org=imget_pixignr( &imhin );
  if (pixignr==(float)INT_MIN)
    {
      pixignr=pixignr_org;
    } 

  npx=imhin.npx;
  npy=imhin.npy;
  
  g=imc_read(&imhin,fpin,0,npx-1,0,npy-1);
  ndat=npx*npy;

  if (pixignr!=pixignr_org)
    {
      for(i=0;i<ndat;i++)
	if(g[i]==pixignr_org) g[i]=pixignr;
    }


  imh_inherit(&imhin,&icomin,&imhout,&icomout,fnamout);

  imhout.npx=npx;
  imhout.npy=npy;
  imhout.ndatx=imhout.npx;
  imhout.ndaty=imhout.npy;
  
  if (bzero==FLT_MAX)
    imc_fits_get_dtype( &imhout, NULL, &bzero, NULL, NULL);
  if (bscale==FLT_MAX)
    imc_fits_get_dtype( &imhout, NULL, NULL, &bscale, NULL);
  if(dtype[0]=='\0'||(fitsdtype=imc_fits_dtype_string(dtype))==-1)
    (void)imc_fits_get_dtype( &imhout, &fitsdtype, NULL, NULL, NULL);
  /* re-set bzero & bscale */
  if(imc_fits_set_dtype(&imhout,fitsdtype,bzero,bscale)==0)
    {
      fprintf(stderr,"Error: Cannot set FITS %s\n",fnamout);
      fprintf(stderr,"Type %d BZERO %f BSCALE %f\n",fitsdtype,bzero,bscale);
      exit(-1);
    }
  imset_pixignr( &imhout, &icomout, pixignr); 

  /*************************************************************/
  if (meshsiz_x==0)
    {
      if (meshsiz_xy>0) meshsiz_x=meshsiz_xy;
      else meshsiz_x=100;
    }
  if (meshsiz_y==0)
    {
      if (meshsiz_xy>0) meshsiz_y=meshsiz_xy;
      else meshsiz_y=100;
    }

  printf("debug: meshsize %d x %d\n",
	 meshsiz_x,meshsiz_y);


  imupdate_fitsf(&icomout, "SSBMSH-I", "%20d / %-46.46s", 
		 meshsiz_x , "SKY MESH X PIXEL" );
  imupdate_fitsf(&icomout, "SSBMSH-J", "%20d / %-46.46s", 
		 meshsiz_y , "SKY MESH Y PIXEL" );
  imaddhistf(&icomout,
	     "Sky subtracted from %s by skysb4a1_spcam",
	     fnamin);

  /* 2008/08/28 */
  for(i=0;i<4;i++)
    {
      sprintf(key,"S_EFMN%1d1",i+1);
      imget_fits_int_value(&icomin,key,&xmin);
      sprintf(key,"S_EFMX%1d1",i+1);
      imget_fits_int_value(&icomin,key,&xmax);
      sprintf(key,"S_EFMN%1d2",i+1);
      imget_fits_int_value(&icomin,key,&ymin);
      sprintf(key,"S_EFMX%1d2",i+1);
      imget_fits_int_value(&icomin,key,&ymax);
      /* copy, not smart at all.. must be rewritten in future*/

      npx2=xmax-xmin+1;
      npy2=ymax-ymin+1;
      h=(float*)malloc(npx2*npy2*sizeof(float));

      printf("debug: region:%d npx2=%d npy2=%d\n",i,npx2,npy2);

      for(iy=ymin-1;iy<ymax;iy++)
	{
	  memcpy(h+(iy-ymin+1)*npx2,g+iy*npx+(xmin-1),npx2*sizeof(float));
	}
      skysub2(h,npx2,npy2,meshsiz_x,meshsiz_y,pixignr,refill,
	      fillf, fillsize_x, fillsize_y,
	      &msky, &msigm
	      );
      for(iy=ymin-1;iy<ymax;iy++)
	{
	  memcpy(g+iy*npx+(xmin-1),h+(iy-ymin+1)*npx2,npx2*sizeof(float));
	}
      free(h);

      imaddhistf(&icomout,
		 " Ch%d  mesh %3dx%3d sky%7.2f+%7.2f",
		 i+1,meshsiz_x,meshsiz_y,msky,msigm);
      sprintf(key,"SKYM%d",i+1);
      imupdate_fitsf(&icomout,key, "%20.4f / %-46.46s", 
		     msky, "TYPICAL SKY ADC ESTIMATED" );
      sprintf(key,"SSGMM%d",i+1);
      imupdate_fitsf(&icomout,key, "%20.4f / %-46.46s", 
		     msigm, "TYPICAL SKY SIGMA ADC ESTIMATED" );
    }
  imclose(fpin,&imhin,&icomin);


  if((fpout=imwopen(fnamout,&imhout,&icomout))==NULL) 
    {
      printf("Cannot open file %s !!\n",fnamout);
      exit(-1);
    }
  imwall_rto( &imhout, fpout, g );
  free(g);
  imclose(fpout,&imhout,&icomout);

  return 0;
}



