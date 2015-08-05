#include <math.h>
#include<stdio.h>


/* */ 
void oi_first(int *no, 
              double *innov,
              double *SRinv,
              double *vec, double *vec1) {
  int nov = no[0];
  int i,j,k;

  k=0;
  for (i=0;i<nov;i++) {
    vec[i]=0;
    vec1[i]=0;
    for (j=0;j<nov;j++) {
/*      if (i==0 && j<10) printf("%9.4f \n",SRinv[k]);*/
      vec[i]=vec[i]+SRinv[k]*innov[j];
      vec1[i]=vec1[i]+SRinv[k];
      k++;
    }
  }
  return;
}
