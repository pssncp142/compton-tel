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
  int n[2] = {3,3};
  double stat[10];
  double event[100];
  int verb = 1;
  double en[2] = {900,1000};
  double vec1[5];
  double vec2[5];
  double cross[3];
  double angle;
  int rnd;

  n_ev = read_bin(argc[1],data,n_proc);
  all_detected(data,n_proc,&n_ev);
  nofproc(data,n_proc,&n_ev,n);
  rnd = floor(((double)rand()/RAND_MAX)*n_ev-1);
  pick_event(event,data,n_proc,n_ev,rnd,verb);
 
  for(i=0;i<event[0];i++){
    for(j=0;j<event[0];j++){
      for(k=0;k<event[0];k++){
	if(k!=j && i!=j && i!=k){
	  for(l=0;l<3;l++) vec1[l] = event[2+4*j+l] - event[2+4*i+l];
	  for(l=0;l<3;l++) vec2[l] = event[2+4*k+l] - event[2+4*i+l];
	  angle = compt_angle(event[1+4*k],event[1+4*k]+event[1+4*i]); 
	  printf("%d %d %d %8.3f %8.3f\n",j,i,k,
		 acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2))),
		 angle);

	}
      }
    }
  }

  free(data);
  free(n_proc);
  return 0;
}
