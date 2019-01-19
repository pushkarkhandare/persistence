#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

void init(double * rx, double * ry) {
    int i = 0;

    FILE *infile;
    infile = fopen("triang.dat", "r");

    while(i<371)
    {
        fscanf(infile,"%lf %lf\n",&rx[i],&ry[i]);
        i++;
    }

    fclose(infile);
}
int main(){
	double * rx;
	double * ry;
	int i = 0;
	int N = 371;
	rx = (double *)malloc(N*sizeof(double));
	ry = (double *)malloc(N*sizeof(double));
	init(rx,ry);
	for (i=0;i<N;i++)
	{
		printf("%lf %lf %d\n",rx[i],ry[i],i);
	}
	printf("%d",time(NULL));
}
