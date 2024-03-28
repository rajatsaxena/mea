//You can include any C libraries that you normally use
#include "stdio.h"
#include "stdlib.h"
#include "omp.h"

#define nthreads 8
#define nbins 101
#define total_length 808

/* 
- written by @Manish - 3/8/17
- We don't need to store all the differences, since we just have to bin it eventually.
  Hence as soon as we find a legitimate difference, we put it in the the corresponding bin.
- As a result we can keep reusing the array ISI[sz]
- The bins go from -1005 to 1005 with a step of 20, so the total number of bins is 101
- The constant "nthreads" is the number of threads, "nbins" is the total number of bins in the histogram
  and "length" is the product of "nbins" and "nthreads"
*/

void binnedISI(double *spikeTimestamps, long *ISICounts, int spikeTimestampsLength) {
  omp_set_num_threads(nthreads);
  int longBin[total_length] = {0};
  #pragma omp parallel
  {
    int i,binIndex;
    int ithread = omp_get_thread_num();
    #pragma omp for
    for (i = 0; i < spikeTimestampsLength; i++) {
      int j;
      for (j = 0; j <= i; j++) {
        //The equality case
        if (i == j) {
          longBin[(ithread*nbins) + 50] += 1;
          //ISICounts[50] += 1;
        }
        else {
          //Calculating the difference between timestamps while converting the difference into milliseconds from microseconds
          double diffTimestamp = (spikeTimestamps[j] - spikeTimestamps[i])/1000.0;
          //remove points greater than 1 sec and less than -1sec
          if (fabs(diffTimestamp) <= 1000.0) {
            //increment the appropriate bins
            binIndex = (int)(diffTimestamp + 1005)/20;
            longBin[(ithread*nbins) + binIndex] += 1;
            longBin[(ithread*nbins) + 100 - binIndex] += 1;
            // ISICounts[binIndex] += 1;
            // ISICounts[100 - binIndex] += 1;
          }
        }
      }
    }
    #pragma omp for
    for (i = 0; i < nbins; i++) {
      int t;
      for (t = 0; t < nthreads; t++) {
        ISICounts[i] += longBin[(nbins*t) + i];
      }
    }
  }
}
  