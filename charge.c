#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "charge.h"

void charge_a(int me,int n,int np ,int *iBeg, int *iEnd)
{
*iBeg=me*(n/np);
*iEnd=(me+1)*(n/np)-1;
}
void charge_b(int me,int n,int np ,int *iBeg, int *iEnd)
{
    int r;
    r=n%np;
    if (me==np-1) {
        *iBeg=me*(n/np);
        *iEnd=n-1;  

    }
    else{
*iBeg=me*(n/np);
*iEnd=(me+1)*(n/np)-1;
    }

}
void charge_c(int me,int n,int np ,int *iBeg, int *iEnd)
{
  int r;
    r=n%np;
    if (me<r) {
*iBeg=me*(n/np+1);
*iEnd=(me+1)*(n/np+1)-1;
    }
    else{
*iBeg=r*(n/np+1)+(me-r)*(n/np);
*iEnd=*iBeg+(n/np)-1; 
    }
}