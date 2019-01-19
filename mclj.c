#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include <time.h>

double total_e ( double * rx, double * ry, int N, double L, double H, double rc, double r0, double sigma, double delta, double gamma) 
{
   int i,j;
   double dx, dy, cx, cy, r, k, r6i;
   double e = 0.0;

   k = r0 - sigma;
   for (i=0;i<(N-1);i++) 
   {
     for (j=i+1;j<N;j++) 
     {
	dx  = (rx[i]-rx[j]);
	dy  = (ry[i]-ry[j]);
	
	cy = round(dy/H);
	cx = dx - cy*gamma*L;
	dx = cx - round(cx/L)*L;
	dy = dy - cy*H;
	
	r = sqrt(dx*dx + dy*dy);
	if (r<sigma) 
	{
	  r6i = (sigma*sigma*sigma*sigma*sigma*sigma)/(r*r*r*r*r*r);
	  e += (r6i*r6i) - (r6i);
	  //printf("%f %f\n",e,r);
	}
	if (r>=sigma && r<=rc)
	{
	  e += (delta*(r - r0)*(r - r0)/(k*k)) - delta;
        }
	//printf("%f %f\n",e,r);
   }
   }
   return e;
}

int main( int argc, char * argv[] ) 
{
  double * rx, * ry;
  int N,c,a,d;
  double L,H;
  double T=1.0,sigma=0.5, r0, delta=1.0, rc,gamma;
  double E_new, E_old;
  double dr=0.1,dx,dy;
  double rxold,ryold;
  int i,j,l,g,s,stepcount=0,m=0,n=0;
  int nCycles = 10000, nSamp, nEq=0;
  int nAcc;
  int short_out=0;
  char *input,buffer1[100],buffer2[100];
  FILE *outptr,*infile,*outconf;
  unsigned long int Seed;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);

  /* Here we parse the command line arguments */
  for (i=1;i<argc;i++) 
  {
    if (!strcmp(argv[i],"-input")) input=argv[++i];
    else if (!strcmp(argv[i],"-sigma")) sigma=atof(argv[++i]);
    else if (!strcmp(argv[i],"-r0")) r0=atof(argv[++i]);
    else if (!strcmp(argv[i],"-delta")) delta=atof(argv[++i]);
    else if (!strcmp(argv[i],"-gamma")) gamma=atof(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dr")) dr=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ne")) nEq = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-so")) short_out=1;
    else if (!strcmp(argv[i],"-seed")) 
      Seed = (unsigned long)atoi(argv[++i]);
    else 
    {
      fprintf(stderr,"Error.  Argument '%s' is not recognized.\n",argv[i]);
      exit(-1);
    }
  }
  infile = fopen(input, "r");
  /*
  stepcount=0;
  while(!feof(infile))
  {
     // fscanf(infile,"%lf %lf",&rx[i],&ry[i]);
      fscanf(infile,"%*[^\n]\n");
      stepcount++;
  }
  rewind(infile);
  N = stepcount;
  */
  fscanf(infile,"%lf %lf %i %lf",&L,&H,&N,&r0);
  printf("%lf %lf %i %lf",L,H,N,r0);
  rx = (double*)malloc(N*sizeof(double));
  ry = (double*)malloc(N*sizeof(double));
  n=0;
  while(!feof(infile))
  {
	  fscanf(infile,"%lf %lf",&rx[n],&ry[n]);
	  n++;
  }
  fclose(infile);

  rc = 2.0*r0 - sigma;
  /* For computational efficiency, use reciprocal T */
  T = 1.0/T;
  nCycles+=nEq;

  for(g=1;g<2;g++)
  {
  gamma = g*0.0001;
  Seed = time(NULL);

  sprintf(buffer1,"/root/persistence/config/lat_%f",gamma);
  sprintf(buffer2,"/root/persistence/energy/energy_%f",gamma);
  outconf = fopen(buffer1,"w");
  outptr = fopen(buffer2,"w");

  /*Shearing*/
  for (s=0;s<N;s++)
  {
	  rx[s] = rx[s] + gamma*ry[s];
  }

  /* Seed the random number generator */
  gsl_rng_set(r,Seed);

  E_old = total_e(rx,ry,N,L,H,rc,r0,sigma,delta,gamma);

  nAcc = 0;
  nSamp = 0;

  for (c=0;c<nCycles;c++)
  {
	 for (d=0;d<N;d++)
	 {
   	  /* Randomly select a particle */
   	  i=(int)gsl_rng_uniform_int(r,N);
   	  /* calculate displacement */
   	  dx = dr*(0.5-gsl_rng_uniform(r));
   	  dy = dr*(0.5-gsl_rng_uniform(r));

   	  /* Save the current position of particle i */
   	  rxold=rx[i];
   	  ryold=ry[i];

   	  /* Displace particle i */
   	  rx[i]+=dx;
   	  ry[i]+=dy;

   	  /* Get the new energy */
   	  E_new = total_e(rx,ry,N,L,H,rc,r0,sigma,delta,gamma);
   	 
	  /* Conditionally accept... */
   	  if (gsl_rng_uniform(r) < exp(-T*(E_new-E_old))) 
	  {
   	    E_old=E_new;
   	    nAcc++;
   	  }
   	  /* ... or reject the move; reassign the old positions */
   	  else 
	  {
   	    rx[i]=rxold;
   	    ry[i]=ryold;
   	  }
         } 
	 fprintf(outptr,"%f\n",E_old);
	 nSamp++;
  }
  fclose(outptr);
 for (m=0;m<N;m++)
 {
	rx[m] = rx[m] - round(rx[m]/L)*L - round(ry[m]/H)*gamma*L;
	ry[m] = ry[m] - round(ry[m]/H)*H;
	fprintf(outconf,"%lf %lf\n",rx[m],ry[m]);
 }
 fclose(outconf);
  /* Output delta-r, the acceptance ratio, 
     and the average energy/particle */
 // if (short_out)
   // fprintf(stdout,"%.6lf %.5lf %.5lf %.5lf\n",
//	    dr,((double)nAcc)/(N*nCycles),
//	    esum/nSamp/N,vir_sum/3.0/nSamp/V+rho*T+pcor);
 // else
 fprintf(stdout,"NVT Metropolis Monte Carlo Simulation"
	    " of the LJ+Harmonic fluid.\n"
	    "---------------------------------------------\n"
	    "Number of particles:              %i\n"
	    "Number of cycles:                 %i\n"
	    "Soft-core radius:                 %.5lf\n"
	    "Equilibrium Bond Length:          %.5lf\n"
	    "Depth of potential:               %.5lf\n"
	    "Cutoff radius:                    %.5lf\n"
	    "Gamma:                            %.5lf\n"
	    "Maximum displacement:             %.5lf\n"
	    "Temperature:                      %.5lf\n"
	    "Results:\n"
	    "Acceptance ratio:                 %.5lf\n"
	    ,
	    N,nCycles,sigma,r0,delta,rc,gamma,dr,1.0/T,((double)nAcc)/(N*nCycles));
}
}
