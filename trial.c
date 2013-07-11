#include "stdio.h"
#include "binio.h"
#include "exclude.h"

int main(){

  double data[50000];
  int n_proc[1000];
  int n_ev;
  int n[2] = {2,3};
  double en[2] = {900,1000};
  char fname[] = "data.bin"; 

  n_ev = read_bin(fname,data,n_proc);
  //clean_0(data,n_proc,&n_ev);
  //nofproc(data,n_proc,&n_ev,n);
  all_detected(data,n_proc,&n_ev);
  //energy_int(data,n_proc,&n_ev,en);
  //proc_exc(data,n_proc,&n_ev,"c");
  print_all(data,n_proc,n_ev);

  return 0;
}
