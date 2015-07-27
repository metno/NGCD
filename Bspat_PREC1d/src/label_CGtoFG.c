#include <math.h>
#include<stdio.h>


void label_CGtoFG(int *nint, double *xlabCG, double *xlabFGtmp, double *xlabFG)
{

  int n = nint[0];
  int i;
  /* cycle over all elements */
  for (i=0; i<n; i++) {
    /* label in FG is missing */
    if (xlabFGtmp[i]==-999.) {
      xlabFG[i]=-999.;
    /* label in CG is defined */
    } else if (xlabCG[i]!=-999.) {
      xlabFG[i]=xlabCG[i];
    /* label in CG is missing and label in FG is defined */
    } else {
      xlabFG[i]=0;
    }
  }
  return;
}
