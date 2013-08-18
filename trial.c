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
#include "reconst.h"

int main(int argv,char **argc){

  srand(time(NULL));
  int i,j,k,l;

  double *data = (double*)malloc(30000000*sizeof(double));
  int *n_proc = (int*)malloc(100000*sizeof(int));
  double *image = (double*)malloc(70000*sizeof(double));
  int *path = (int*)malloc(70000000*sizeof(int));  
  double *cones = (double*)malloc(1000000*sizeof(double));
  for(i=0;i<70000;i++) image[i]=0;
  double cone[15];
  //double cones[10000000]={0};
  int n_ev;
  int n[2] = {3,7};
  double stat[10];
  double event[200];
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
    //printf("%d\n",compt_analyse(cone,event));
    if(compt_analyse(cone,event)==1 && cone[12] < 0.01){
    //printf("%d\n",cs++);
      rec_backproj(image,cone,5000,0.1*3.14,100);
      //printf("b\n");
    }
    //printf("%f\n",event[0]);
    /*if(event[0] > 2 && event[0] < 7){
      //printf("%f\n",event[0]);
      st = compt_get_path(path,event,verb);
      if(st == 0){
	//back_add_cone(cones,event,path,1,verb);
	back_add_cone_1(cone,event,path);
	rec_backproj(image,cone,5000,0.1*3.14,100);
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

  //back_proj(cones,5000,0.1*3.14,100,verb);
  
  //free(data);
  //free(n_proc);
  
  FILE* f = fopen("image.txt","w+");
  printf("%f \n",image[50000]);
  
  for(i=0;i<100;i++){
    for(j=0;j<100;j++){
      //printf("%d %d %d\n",i,j,i*256+j);
      fprintf(f,"%f ",image[j*100+i]);
    }
    fprintf(f,"\n");
  }
  fclose(f);

  //free(image);
  return 0;
}
