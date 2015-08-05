#include <math.h>
#include<stdio.h>


/* bilinear interpolation from a coarser grid (cg) to a finer grid (fg) */
void oi_xb_fast(int *ng, int *no, 
                double *xg, double *yg, double *zg,
                double *xo, double *yo, double *zo,
                double *Dh, double *Dz,
                double *vec1,
                double *bparam,
                double *xb, 
                double *xidi) {

  int ngv = ng[0];
  int nov = no[0];
  double Dhv = Dh[0];
  double Dzv = Dz[0];
  int i,j;
  double g,Dh2,Dz2,hd2,vd2;
  double zinv,zabov,zbelo,xb_ab,xb_bl;
  double h0,h1;

  Dh2=Dhv*Dhv;
  Dz2=Dzv*Dzv;
  
  for (i=0;i<ngv;i++) {
    if (bparam[11]==0) {
      xb[i]=bparam[2]+bparam[5]*(xg[i]-bparam[9])+
                       bparam[7]*(yg[i]-bparam[10])+
                       bparam[3]*zg[i];
    } else if (bparam[11]==1) {
      zinv=bparam[0];
      zabov=zinv+bparam[1];
      zbelo=zinv-bparam[1];
      if (zg[i]>zabov) {
        xb[i]=bparam[2]+bparam[5]*(xg[i]-bparam[9])+
                         bparam[7]*(yg[i]-bparam[10])+
                         bparam[3]*(zg[i]-zinv);
      } else if (zg[i]<=zbelo) {
        xb[i]=bparam[2]+bparam[6]*(xg[i]-bparam[9])+
                         bparam[8]*(yg[i]-bparam[10])+
                         bparam[4]*(zg[i]-zinv);
      } else {
        xb_ab=bparam[2]+bparam[5]*(xg[i]-bparam[9])+
                         bparam[7]*(yg[i]-bparam[10])+
                         bparam[3]*(zg[i]-zinv);
        xb_bl=bparam[2]+bparam[6]*(xg[i]-bparam[9])+
                         bparam[8]*(yg[i]-bparam[10])+
                         bparam[4]*(zg[i]-zinv);
        xb[i]=(xb_ab*(zg[i]-zbelo)+xb_bl*(zabov-zg[i]))/(zabov-zbelo);
      }
    } else if (bparam[11]==2) {
      h0=bparam[0];
      h1=bparam[0]+bparam[1];
      if (zg[i]>=h1) {
        xb[i]=bparam[2]+bparam[3]*zg[i]+
                         bparam[5]*(xg[i]-bparam[9])+
                         bparam[7]*(yg[i]-bparam[10]);
      } else if (zg[i]<=h0) {
        xb[i]=bparam[2]+bparam[3]*zg[i]-bparam[4]+
                          bparam[6]*(xg[i]-bparam[9])+
                          bparam[8]*(yg[i]-bparam[10]);
      } else {
        xb[i]=bparam[2]+bparam[3]*zg[i]-bparam[4]/2.*(1+cos(M_PI*(zg[i]-h0)/bparam[1])) + 
               ( (bparam[5]*(xg[i]-bparam[9])+bparam[7]*(yg[i]-bparam[10]))*(zg[i]-h0) + 
                 (bparam[6]*(xg[i]-bparam[9])+bparam[8]*(yg[i]-bparam[10]))*(h1-zg[i]) ) / bparam[1];
      }
    }
    xidi[i]=0;
  }

  for (j=0;j<nov;j++) { 
    for (i=0;i<ngv;i++) {
      hd2=((xg[i]-xo[j])*(xg[i]-xo[j])+(yg[i]-yo[j])*(yg[i]-yo[j])) / (1000.*1000.);
      vd2=(zg[i]-zo[j])*(zg[i]-zo[j]);
      g=exp(-0.5*(hd2/Dh2+vd2/Dz2));
      xidi[i]=xidi[i]+g*vec1[j];
    }
  }

  return;
}
