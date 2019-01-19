#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double autocorr(double * somedata, int t, int tmax)
{
	int i,j;
	double corr=0;

	for (i=0;i<(tmax - t);i++)
	{
		corr += somedata[i+t]*somedata[i];
	}
	corr = corr/tmax;
	return corr;
}
int main(int argc, char * argv[])
{
	char *input,buffer[100];
	double corr1,corr2,corr3;
	double * somedata;
	int stepcount,n,t,tmax;
	FILE *infile;
	
	input = argv[1];
	infile = fopen(input, "r");
  	stepcount=0;
  	while(!feof(infile))
  	{
      		fscanf(infile,"%*[^\n]\n");
      		stepcount++;
  	}
  	rewind(infile);
  	tmax = stepcount;
  	somedata = (double*)malloc(tmax*sizeof(double));
  	n=0;
  	while(!feof(infile))
  	{
	  	fscanf(infile,"%lf",&somedata[n]);
	  	n++;
  	}
  	fclose(infile);
	for (t=0;t<tmax;t++)
	{
		

}

