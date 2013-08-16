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
#include "backproj.h"
#include "analyse.h"

int main(int argv,char **argc){

  srand(time(NULL));
  int i,j,k,l;

  double *data = (double*)malloc(300000000*sizeof(double));
  int *n_proc = (int*)malloc(1000000*sizeof(int));
  int *path = (int*)malloc(100*sizeof(int));
  double cone[15];
  int n_ev;
  int n[2] = {3,15};
  double stat[10];
  double event[200];
  double cones[1000000]={0};
  int verb = 0;
  double en[2] = {900,1000};
  double vec1[5];
  double vec2[5];
  double cross[3];
  double angle;
  int rnd;
  int st;
  int cs=0,cf=0,cf2=0,us=0,uf=0;

  n_ev = read_bin(argc[1],data,n_proc);
  //nofproc(data,n_proc,&n_ev,n);
  //clean_out(data,n_proc,&n_ev);

  for(i=0;i<n_ev;i++){
    pick_event(event,data,n_proc,n_ev,i,verb);
    compt_analyse(event,path);
    //printf("%f\n",event[0]);
    /*if(event[0] > 2 && event[0] < 4){
      //printf("%f\n",event[0]);
      st = compt_get_path(path,event,verb);
      if(st == 0){
	back_add_cone(cones,event,path,1,verb);
	cs++;
      } else if(st == 1){
	cf++;
      } else if(st == 2){
	us++;
      } else if(st == 3){
	uf++;
      } else if(st == 4){
	cf2++;
      }
    
      if(verb)
	printf("\ncs :%4d   cf :%4d   cf2 :%4d   us :%4d   uf :%4d\n",
	       cs,cf,cf2,us,uf);
	       }*/
  }

  //back_proj(cones,500,0.2*3.14,100,verb);
  
  free(data);
  free(n_proc);
  //  free(path);
  return 0;
}
