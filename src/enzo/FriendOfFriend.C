/***********************************************************************
/
/  FOF ROUTINE: MERGE PARTICLE USE FRIEND-OF-FRIEND ALGORITHM
/
/  written by: Peng Wang
/  date:       Januaray, 2009
/  modified1: Nikhil Bisht, March 2023
/             rewrote the fof algorithm completely        
/
/  PURPOSE:
/
/  DESCRIPTION:
/
************************************************************************/

#include <math.h>
#include "macros_and_parameters.h"

int fof(FLOAT *x, FLOAT *y, FLOAT *z, const int &np, const FLOAT &l,
        int *group, int &ng)
{
  if (np <= 1)
    return 1;

  for (int i = 0; i < np; i++) {
    if (group[i] == -1){
      group[i] = ng;
      FLOAT l2 = l*l;
      for (int ip = 0; ip < np; ip++) {
        if (ip != i && group[ip] == -1){
          FLOAT r2 = pow(x[ip] - x[i],2) + pow(y[ip] - y[i], 2) + pow(z[ip] - z[i], 2);
          if (r2 <= l2) {
            group[ip] = group[i];
          }
        }
      }
      ng++;
    }
  }
  return 1;
}
