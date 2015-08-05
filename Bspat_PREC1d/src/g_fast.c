#include <math.h>
#include<stdio.h>


/* bilinear interpolation from a coarser grid (cg) to a finer grid (fg) */
void g_fast(int *ng, int *no, 
            double *xg, double *yg, double *zg,
            double *xo, double *yo, double *zo,
            double *Dh, double *Dz,
            double *g) {
  int ngv = ng[0];
  int nov = no[0];
  int i,j,k;
  double Dh2,Dz2,hd2,vd2;

  k=0;
  Dh2=Dh[0]*Dh[0];
  Dz2=Dz[0]*Dz[0];
/*  for (i=0;i<ngv;i++) {
    for (j=0;j<nov;j++) { */
  for (j=0;j<nov;j++) { 
    for (i=0;i<ngv;i++) {
      hd2=((xg[i]-xo[j])*(xg[i]-xo[j])+(yg[i]-yo[j])*(yg[i]-yo[j])) / (1000.*1000.);
      vd2=(zg[i]-zo[j])*(zg[i]-zo[j]);
      g[k]=exp(-0.5*(hd2/Dh2+vd2/Dz2));
      k++;
    }
  }
  return;
}
