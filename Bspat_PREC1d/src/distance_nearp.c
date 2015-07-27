#include <math.h>
#include<stdio.h>


void foo(int *ngi, int *npi, double *xg, double *yg,double *xp,double *yp, double *dist)
{

  int ng = ngi[0];
  int np = npi[0];
  int i,j;
  double aux,dx2,dy2;

  printf("%i\n",ng);
  printf("%i\n",np);

  for (i=0; i<ng; i++) {
    dist[i] = pow( ((xg[i] - xp[0]) * (xg[i] - xp[0]) + (yg[i] - yp[0]) * (yg[i] - yp[0])), 0.5);
    if ((i % 10000) ==0) printf("%i\n",i);
    for (j=1; j<np; j++) {
      aux = pow( ((xg[i] - xp[j]) * (xg[i] - xp[j]) + (yg[i] - yp[j]) * (yg[i] - yp[j])), 0.5);
      if ( aux<dist[i] ) {
        dist[i] = aux;
      }
    }
  }
  return;
}
