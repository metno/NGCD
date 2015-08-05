#include <math.h>
#include<stdio.h>


/* bilinear interpolation from a coarser grid (cg) to a finer grid (fg) */
void cg2fg_cgweights(int *nfg, int *ncg, 
                     double *fg_weight, double *cg_weight,
                     double *cg_mask,
                     int *cg_c1, int *cg_c2,int *cg_c3, int *cg_c4,
                     double *na)
/*Note: the output fg field fg_field[1,..,nfg] is computed only for points where 
  fg_mask[1,..,nfg] is not "na" (if fg_mask[i]==na -> fg_field[i]=na)
description: given that: for each fg point the surrounding 4 cg points 
  (we assume they always exist!) are in:
  cg_c1 (nearest), ..., cg_c4 (far)
 the fg output field is obtained resampling the cg field (cg_field)
 by using the bilinear interpolation.
Note:
cg_wX <- weight to be used for the spatial interpolation (in case of not nas)

example:
+cg1[i]---cg2[i]+
| fg[i]+        | 
|               |
|               |
+cg3[i]---cg4[i]+
if fg_mask[i]!=na then fg[i]=cg_w1[i]*cg1[i]+...+cg_w4[i]*cg4[i]
NOTE: NA values are not allowed in cg!
*/
{
  int nfgv = nfg[0];
  int ncgv = ncg[0];
  double nav = na[0];
  int cg_count[ncgv];
  int i;

  for (i=0; i<ncgv; i++) {
    if (cg_mask[i]!=nav) {
      cg_count[i]=0;
      cg_weight[i]=0;
    } else {
      cg_count[i]=nav;
      cg_weight[i]=nav;
    }
  }

  for (i=0; i<nfgv; i++) {
    if (cg_mask[cg_c1[i]]!=nav) {
      cg_count[cg_c1[i]]=cg_count[cg_c1[i]]+1;
      cg_weight[cg_c1[i]]=cg_weight[cg_c1[i]]+(fg_weight[i]-cg_weight[cg_c1[i]])/cg_count[cg_c1[i]];
    }
  }
  for (i=0; i<ncgv; i++) {
    if (cg_mask[i]!=nav && cg_count[i]==0) cg_weight[i]=1;
  }
  return;
}
