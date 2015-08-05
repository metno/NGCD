#include <math.h>
#include<stdio.h>

/* nearest neighbour interpolation from a coarser grid (cg) to a finer grid (fg) */
void cg2fg_ngb(int *nfg, double *fg_mask, double *fg_field, 
               double *cg_field,
               int *cg_c1, int *cg_c2,int *cg_c3, int *cg_c4,
               double *na, int *na_rm)
/* Note: the output fg field fg_field[1,..,nfg] is computed only for points where 
 fg_mask[1,..,nfg] is not "na" (if fg_mask[i]==na -> fg_field[i]=na)
description: 
 given that: for each fg point the surrounding 4 cg points 
  (we assume they always exist!) are in:
  cg_c1 (nearest), ..., cg_c4 (far)
 the fg output field is obtained resampling the cg field (cg_field) by using the
 nearest neighbour interpolation.
Example:
+cg1[i]---cg2[i]+
| fg[i]+        | 
|               |
|               |
+cg3[i]---cg4[i]+
if (na_rm!=1, i.e. false)
  --> fg[i]=cg1[i]
if (na_rm==0, i.e. true)
  --> fg[i]=cg1[i], if cg1[i] is not na
      fg[i]=cg2[i], if cg2[i] is not na 
      fg[i]=cg3[i], if cg3[i] is not na 
      fg[i]=cg4[i], if cg4[i] is not na 
      fg[i]=na, otherwise
*/
{

  int nfgv = nfg[0];
  double nav = na[0];
  int na_rmv = na_rm[0];
  int i,j,k;

  for (i=0; i<nfgv; i++) {
    fg_field[i]=nav;
    if (fg_mask[i]!=nav) {
      for (j=1; j<=4; j++) {
        if (j==1) fg_field[i]=cg_field[cg_c1[i]];
        if (j==2) fg_field[i]=cg_field[cg_c2[i]];
        if (j==3) fg_field[i]=cg_field[cg_c3[i]];
        if (j==4) fg_field[i]=cg_field[cg_c4[i]];
        if (na_rmv!=1) {
          break;
        } else if (fg_field[i]!=nav) {
          break;
        }
      }
    }
  }
  return;
}
