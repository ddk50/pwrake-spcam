#! /usr/bin/perl

if (@ARGV!=2)
{
    printf("\nUsage: getkey.pl KEYWORD filename\n");
    printf("         eg) getkey.pl NAXIS1 test.fits\n\n");
    exit(-1);
}

$key=$ARGV[0];
$FITS=$ARGV[1];
open(FITS);

while(<FITS>)
{
 $LINE=$_; 
 if( $LINE=~/$key\s*=\s+\'([^\']+[^ \']+) *\'/ || 
     $LINE=~/$key\s*=\s+([0-9\.Ee\+\-]+)\s/ ||
     $LINE=~/$key\s*=\s+([TF])\s/)
 { 	
     printf("%s\n",$1); last;
 }
 elsif ($LINE=~/END {77}/)
 {
     last;
 }
}
