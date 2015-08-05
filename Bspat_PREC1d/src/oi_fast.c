#include <math.h>
#include<stdio.h>


/* bilinear interpolation from a coarser grid (cg) to a finer grid (fg) */
void oi_fast(int *ng, int *no, 
             double *xg, double *yg, double *zg,
             double *xo, double *yo, double *zo,
             double *Dh, double *Dz,
             double *xb, 
             double *vec,
             double *vec1,
             double *smooth,
             double *xinf,
             double *xa, double *xidi) {
  int ngv = ng[0];
  int nov = no[0];
  double xinfv = xinf[0];
  double Dhv = Dh[0];
  double Dzv = Dz[0];
  int i,j;
  double g,Dh2,Dz2,hd2,vd2;

  Dh2=Dhv*Dhv;
  Dz2=Dzv*Dzv;
  for (i=0;i<ngv;i++) {
    xa[i]=xb[i];
    xidi[i]=0;
  }
  for (j=0;j<nov;j++) { 
    for (i=0;i<ngv;i++) {
      hd2=((xg[i]-xo[j])*(xg[i]-xo[j])+(yg[i]-yo[j])*(yg[i]-yo[j])) / (1000.*1000.);
      vd2=(zg[i]-zo[j])*(zg[i]-zo[j]);
      g=exp(-0.5*(hd2/Dh2+vd2/Dz2));
      xa[i]=xa[i]+g*vec[j];
      xidi[i]=xidi[i]+g*vec1[j];
    }
  }
  for (i=0;i<ngv;i++) {
    xa[i]=xa[i]*smooth[i];
    if (xa[i]<xinfv && smooth[i]>=0.9) xa[i]=xinfv; 
    if (xa[i]<0) xa[i]=0.;
  }

  return;
}
