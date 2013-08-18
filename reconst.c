#include "math.h"
#include "3d_cart_vec.h"
#include "stdlib.h"
#include "time.h"
#include "stdio.h"

#define PI 3.14f

double rec_angle(double post, double pre){
  double rest_e=511;
  return acos(1-rest_e*(1/post-1/pre));
}

int rec_backproj(double *image, double *cone, double height, double ang, int n_grid){

  int i,j,k;
  double posx[70000];
  double posy[70000];
  double length = height*tan(ang);
  double g_size = length*2/n_grid;
  double g_size_over_2 = g_size*0.5;
  double c_z[3],c_y[3],c_x[3],dir[3],d_x[3]={1,0,0};
  double phi;
  double ratio;
  double theta;

  for(i=0;i<n_grid;i++){
    for(j=0;j<n_grid;j++){
      posx[j*n_grid+i] = (i+0.5)*g_size-length;
      posy[j*n_grid+i] = (j+0.5)*g_size-length;
    }
  }

  c_z[0] = cone[1]-cone[4];
  c_z[1] = cone[2]-cone[5];
  c_z[2] = cone[3]-cone[6];
  vec_unit(c_z);
  theta = rec_angle(cone[11],cone[10]);
  if(theta>PI*0.5) {
    theta = PI-theta;
    vec_scap(-1,c_z);
  }
  vec_ort_wc(c_x,d_x,c_z);
  vec_cross_wc(c_y,c_x,c_z);
  vec_rotate(c_y,PI*0.5-theta,c_x);
  for(j=0;j<5000;j++){
    phi = (double) rand()/RAND_MAX;
    phi *= 2*PI;
    vec_rotate_wc(dir,c_z,phi,c_x);
    ratio = height/dir[2];
    dir[0] *= ratio;
    dir[1] *= ratio;
    if(dir[2] > 0 && dir[0] < length && dir[0] > -length && dir[1] < length && dir[1] > -length){
      //image[(int)floor((dir[0]+length)/g_size-0.5)*n_grid+(int)floor((dir[1]+length)/g_size-0.5)]=cone[10]/5000.;
      for(k=0;k<n_grid*n_grid;k++){
	if(dir[0] <= posx[k]+g_size_over_2 && dir[0] > posx[k]-g_size_over_2 &&
	   dir[1] <= posy[k]+g_size_over_2 && dir[1] > posy[k]-g_size_over_2){
	  image[k] += cone[10]/5000.; break;
	}
      }
    }
  }

  return 0;

}
