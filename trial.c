#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "binio.h"
#include "exclude.h"
#include "event.h"
#include "stat.h"
#include "3d_cart_vec.h"
#include "compton.h"

int main(int argv,char **argc){

  srand(time(NULL));
  int i,j,k,l;
  double *data = (double*)malloc(1000000*sizeof(double));
  int *n_proc = (int*)malloc(10000*sizeof(int));
  int n_ev;
  int n[2] = {4,4};
  double stat[10];
  double event[100];
  double match[1000];
  double path[1000];
  int verb = 0;
  double en[2] = {900,1000};
  double vec1[5];
  double vec2[5];
  double cross[3];
  double angle;
  int rnd;

  n_ev = read_bin(argc[1],data,n_proc);
  all_detected(data,n_proc,&n_ev);
  nofproc(data,n_proc,&n_ev,n);
  for(i=0;i<n_ev;i++){
    pick_event(event,data,n_proc,n_ev,i,verb);
    compt_match3(match,path,event,verb);
  }

  free(data);
  free(n_proc);
  return 0;
}
