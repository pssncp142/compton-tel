#include "stdio.h"
#include "stdlib.h"
#include "binio.h"
#include "exclude.h"
#include "stat.h"

int main(){

  double *data = (double*)malloc(1000000*sizeof(double));
  int *n_proc = (int*)malloc(10000*sizeof(int));
  int n_ev;
  int n[2] = {8,8};
  double stat[10];
  double en[2] = {900,1000};
  char fname[] = "data/800-10000.bin"; 

  n_ev = read_bin(fname,data,n_proc);
  //proc_exc(data,n_proc,&n_ev,"p");
  nofproc(data,n_proc,&n_ev,n);
  stat_complete(stat,data,n_proc,n_ev,1);
  //print_all(data,n_proc,n_ev);

  free(data);
  free(n_proc);
  return 0;
}
