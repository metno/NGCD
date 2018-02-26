#include <math.h>
#include<stdio.h>


/* optimal interpolation */
void oi_fast(int *ng,
             int *no, 
             double *xg, 
             double *yg, 
             double *zg,
             double *lg,
             double *xo, 
             double *yo, 
             double *zo,
             double *lo,
             double *Dh, 
             double *Dz,
             double *lafmin,
             double *xb, 
             double *vec,
             double *vec1,
             double *xa, 
             double *xidi) {
  int ngv = ng[0];
  int nov = no[0];
  double Dhv = Dh[0];
  double Dzv = Dz[0];
  double lafmin0 = lafmin[0];
  int i,j;
  double g,Dh2,Dz2,hd2,vd2,ld;

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
      ld=1-(1-lafmin0)*(lg[i]-lo[j])*(lg[i]-lo[j]);
      g=exp(-0.5*(hd2/Dh2+vd2/Dz2))*ld;
      xa[i]=xa[i]+g*vec[j];
      xidi[i]=xidi[i]+g*vec1[j];
    }
  }

  return;
}
