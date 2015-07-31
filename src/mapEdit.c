#include "R.h"
/* several routines expect integer values for x/y! */

void mapvalence(int*x,int*y,int*len,int*valence);
void mapclean(int*,int*,int*,int*,int*,int*);
void mapsplit(int*x,int*y,int*len,
               int*nx,int*ny,int*nlen,
               int*valence,int*gon,int*ngon);

void mapclean(int*x,int*y,int*len,int*nx,int*ny,int*nlen){
  int i,j;
/* BUGfix: mini-islands that are 1 point repeated twice become corrupted! */ 
/* if there is only one point in the polygon, you want it repeated! */
  nx[0]=x[0];
  ny[0]=y[0];
  j=1;
  for (i=1 ; i<*len ; i++){
    if( ISNA(x[i]) || ISNA(nx[j-1]) ||
        (j>1 && ISNA(nx[j-2]) && ISNA(x[i+1])) ||
        x[i]!=nx[j-1] || y[i]!=ny[j-1]){
       nx[j]=x[i];
       ny[j]=y[i];
       j++;
    }
  }
  *nlen=j;
}

void mapvalence(int*x,int*y,int*len,int*valence){
  int i,j,val;
  
  for(i=0 ; i<*len ; i++){
    if(x[i] == NA_INTEGER) valence[i] = 0;
    else{
       val=0;
       for(j=0 ; j<i ; j++) {
         if( (x[j]!=NA_INTEGER) && !(x[i]-x[j]) && !(y[i]-y[j]) ) val = valence[j];
       }
       if(!val) {
         for(j=i; j<*len;j++)  
         if( (x[j]!=NA_INTEGER) && !(x[i]-x[j]) && !(y[i]-y[j]) ) val++;
       }
       valence[i] = val;
    }
  }
}

void mapsplit(int*x,int*y,int*len,
               int*nx,int*ny,int*nlen,
               int*valence,int*gon,int*ngon){
/* BUG: coastlines of length 1 (both points are vertices) are not identified */
  int i,ip,igon,ilin;

  nx[0]=x[0];
  ny[0]=y[0];
  ip=1;
  ilin=1;
  igon=0;
  for(i=1;i<*len;i++){
    nx[ip]=x[i];ny[ip]=y[i];
    ip++;
    if(valence[i]==0) { /* end/start of a polygon */
      gon[igon]=ilin;
      igon++;ilin++;
      gon[igon]=NA_REAL;
      igon++;
    }
    if( (valence[i]==2 && valence[i+1]==1 && valence[i-1]>0) || 
        (valence[i]==2 && valence[i-1]==1 && valence[i+1]>0) ||
        (valence[i] > 2 && valence[i-1]>0 && valence[i+1]>0 )) {
      nx[ip]=ny[ip]=NA_REAL;
      ip++;
      nx[ip]=x[i];ny[ip]=y[i];
      ip++;
      gon[igon]=ilin;
      igon++;ilin++;
    }
  }

  *nlen=ip;
  *ngon=igon;
}


