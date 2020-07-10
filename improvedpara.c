#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "utils.h"

double work_it_par(long *old, long *new, long *super, long *simple, long *fibonacci) {
  int i, j, k;
  int u, v, w;
  int ton = 0;
  long compute_it, moving_average;
  double pi, pi2, x , y, sum, step = 0.0;
  long dot_product=0;
  long local_nCirc, nCirc=0;
  long aggregate=1.0;
  double r=1.0;
  int was_smart = 16;
  int threadid;
  long DIMSQ = DIM * DIM;
  long km;
  long kp; 
  long idimsq, ipdimsq;
  long imdimsq;
  long jdim, jmdim, jpdim;
  long temp, temp2;
  long tempvar;

//#pragma omp parallel
//#pragma omp set_num_threads(2)
 {
      //#pragma omp for schedule(dynamic, 1)
      for(i=0; i<DIM-1;i++)
      {
         super[i] += simple[i];
      }
 }

//#pragma omp parallel shared(dot_product) private(moving_average)
  {
  // #pragma omp for schedule(static,1) reduction(+:dot_product)
  // #pragma omp for reduction(+:dot_product)   
   for(i=0; i<DIM-1;i++)
      {
        dot_product += super[i]*simple[i];
       
        moving_average = 0; 
        for(ton=i;ton<DIM-1-WINDOW_SIZE;ton++)
        {
          moving_average += simple[ton];
        }
      }
  }

  int a_secret = 5;
  fibonacci[0] = 1;
  fibonacci[1] = 1;

  for(i=2; i<DIM-1;i++)
  {
    fibonacci[i]=fibonacci[i-1]+fibonacci[i-2];
    if(i==3)
    {
       printf("\n A secret is: %d",obfuscate_obfuscate_obfuscate(a_secret));
    }
  }


  step = 1.0 / NUM_STEPS; 
  for (i=0;i<NUM_STEPS; i++)
  {
    x = (i+0.5)*step;
    sum = sum + 4.0/(1.0+x*x);
  }
  pi = step * sum;

 printf("\n %d trials, Riemann flavored pi is %f \n",NUM_STEPS, pi); 
 
// #pragma omp parallel private(x, y, r, threadid) shared(nCirc)

//#pragma omp parallel private(x,y,r) shared(nCirc)
 {
      // printf("Currently in thread #%d\n",omp_get_thread_num());
      // #pragma omp for reduction(+:nCirc) 
      // #pragma omp for schedule(dynamic, 1) reduction(+:nCirc) 
      for(i=0;i<NUM_TRIALS; i++)
      {
        x = (random()%10000000)/10000000.0; 
        y = (random()%10000000)/10000000.0;
        if (( x*x + y*y) <= r) {
          nCirc++;
        }
      }
 }
  pi2 = 4.0 * ((double)nCirc/(double)NUM_TRIALS);
 
  printf("\n %d trials, Monte-Carlo flavored pi is %f \n",NUM_TRIALS, pi2); 

  long wntf = we_need_the_func();
  long gtf = gimmie_the_func();


#pragma omp parallel 
{
  #pragma omp for private(k, compute_it) reduction(+:aggregate) schedule(dynamic, 1)
  for (i=1; i<DIM-1; i++) {
    long idimsq = i * DIMSQ;
      for (j=1; j<DIM-1; j++) {
      long jdim = j * DIM;
       
         for (k=1; k<DIM-1; k++) {
         compute_it = old[idimsq+jdim+k] * wntf; 
         aggregate+= compute_it / gtf; 
      }
    }
  }
}


printf("AGGR:%ld\n",aggregate);



  for (i=1; i<DIM-1; i++) {
    idimsq = i * DIMSQ;  
    imdimsq = (i-1) * DIMSQ;  
    ipdimsq = (i+1) * DIMSQ;  
    for (j=1; j<DIM-1; j++) {
        jdim = j * DIM;  
        jmdim = (j-1) * DIM;  
        jpdim = (j+1) * DIM;  
      for (k=1; k<DIM-1; k++) {
	   kp = k+1;
	   km = k-1; 
           temp = idimsq+jdim+k;
	   tempvar = 0;
           new[temp]=0;
   	   tempvar =  
		   (old[imdimsq+jmdim+km]
		   +old[imdimsq+jmdim+k]
		   +old[imdimsq+jmdim+kp]

		   +old[imdimsq+jdim+km]
		   +old[imdimsq+jdim+k]
		   +old[imdimsq+jdim+kp]

		   +old[imdimsq+jpdim+km]
		   +old[imdimsq+jpdim+k]
		   +old[imdimsq+jpdim+kp]

		   +old[idimsq+jmdim+km]
		   +old[idimsq+jmdim+k]
		   +old[idimsq+jmdim+kp]
		   
		   +old[idimsq+jdim+km]
		   +old[idimsq+jdim+k]
		   +old[idimsq+jdim+kp]

		   +old[idimsq+jpdim+km]
		   +old[idimsq+jpdim+k]
		   +old[idimsq+jpdim+kp]

		   +old[ipdimsq+jmdim+km]
		   +old[ipdimsq+jmdim+k]
		   +old[ipdimsq+jmdim+kp]

		   +old[ipdimsq+jdim+km]
		   +old[ipdimsq+jdim+k]
		   +old[ipdimsq+jdim+kp]

		   +old[ipdimsq+jpdim+km]
		   +old[ipdimsq+jpdim+k]
		   +old[ipdimsq+jpdim+kp])/27;

	   new[temp]= tempvar;
          
/*


           new[temp]+=old[imdimsq+jmdim+km];
           new[temp]+=old[imdimsq+jmdim+k];
           new[temp]+=old[imdimsq+jmdim+kp];

           new[temp]+=old[imdimsq+jdim+km];
           new[temp]+=old[imdimsq+jdim+k];
           new[temp]+=old[imdimsq+jdim+kp];

           new[temp]+=old[imdimsq+jpdim+km];
           new[temp]+=old[imdimsq+jpdim+k];
           new[temp]+=old[imdimsq+jpdim+kp];

           new[temp]+=old[idimsq+jmdim+km];
           new[temp]+=old[idimsq+jmdim+k];
           new[temp]+=old[idimsq+jmdim+kp];
           
           new[temp]+=old[idimsq+jdim+km];
           new[temp]+=old[idimsq+jdim+k];
           new[temp]+=old[idimsq+jdim+kp];

           new[temp]+=old[idimsq+jpdim+km];
           new[temp]+=old[idimsq+jpdim+k];
           new[temp]+=old[idimsq+jpdim+kp];

           new[temp]+=old[ipdimsq+jmdim+km];
           new[temp]+=old[ipdimsq+jmdim+k];
           new[temp]+=old[ipdimsq+jmdim+kp];

           new[temp]+=old[ipdimsq+jdim+km];
           new[temp]+=old[ipdimsq+jdim+k];
           new[temp]+=old[ipdimsq+jdim+kp];

           new[temp]+=old[ipdimsq+jpdim+km];
           new[temp]+=old[ipdimsq+jpdim+k];
           new[temp]+=old[ipdimsq+jpdim+kp];


           new[temp]/=27;
*/

      }
    }
  }

    for (i=1; i<DIM-1; i++) {
    idimsq = i*DIM*DIM;
    for (j=1; j<DIM-1; j++) {
    jdim = j*DIM; 
    for (k=1; k<DIM-1; k++) {
        u=(new[idimsq+jdim+k]/100);      
        if (u<=0) u=0;
        if (u>=9) u=9;
        histogrammy[u]++;
      }
    }
  }

return (double) (dot_product+moving_average+pi+pi2);


}
