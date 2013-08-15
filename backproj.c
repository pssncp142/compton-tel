/*****************************************************************************\
 * Yigit Dallilar 12.07.2013                                                 *
 * DTU-Space                                                                 *
 * Compton telescope library for handling back projection of the cones       *
\*****************************************************************************/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "3d_cart_vec.h"
#include "compton.h"

#define PI 3.14159f

static int i,j,k,l;

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int back_add_cone(double cones[], double event[], int path[], 
		  int time, int verb){
  
  int l_p=path[0];
  int n_c=(int)cones[0];
  double en=0;
  double vec1[3];
  double axisz[3] = {0,0,1}; 
  double axisx[3] = {1,0,0};
  
  if(time>l_p) {
    printf("!!!ERROR : time : %d  and l_p : %d varibles do not check..\n",
	   time,l_p);
    return 1;
  }

  if(verb) printf("\n%d cones will be added :\n",time);

  for(j=0;j<time;j++){

    for(i=0;i<l_p-1;i++){
      en += event[1+4*path[2+i]];
    }
    for(i=0;i<3;i++){
      vec1[i] = event[2+4*path[1+l_p]+i] - event[2+4*path[1+l_p-1]+i];
    }
    cones[1+4*n_c] = 
      acos(vec_dotp(vec1,axisz)/(vec_norm(vec1)*vec_norm(axisz)));
    cones[1+4*n_c+1] = vec_angle(axisz,vec1,axisx);
    cones[1+4*n_c+2] = en + event[1+4*path[1+l_p]];
    cones[1+4*n_c+3] = compt_angle(en,cones[1+4*n_c+2]);

    if(cones[1+4*n_c+3]>PI*0.5){
      cones[1+4*n_c]   = PI - cones[1+4*n_c];
      cones[1+4*n_c+3] = PI - cones[1+4*n_c+3];
    }

    if(verb){
      printf("-Cone Index --> %3d\n",(int)cones[0]);
      printf(" * Cone Axis Theta : %8.3f PI\n",cones[1+4*n_c]/PI);
      printf(" * Cone Axis Phi   : %8.3f PI\n",cones[1+4*n_c+1]/PI);
      printf(" * Cone Top Theta  : %8.3f PI\n",cones[1+4*n_c+3]/PI);
      printf(" * Photon energy   : %8.3f keV\n",cones[1+4*n_c+2]);
    }
    cones[0]++;

  }

  return 0;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int back_proj(double cones[], double height, double ang, 
	      int n_grid, int verb){

  double image[20000]={0};
  double posx[20000];
  double posy[20000];
  double length = height*tan(ang);
  double g_size = length*2/n_grid;
  double g_size_over_2 = g_size*0.5;
  double c_z[3],c_y[3],c_x[3],dir[3],d_x[3]={1,0,0};
  double phi;
  double ratio;

  for(i=0;i<n_grid;i++){
    for(j=0;j<n_grid;j++){
      posx[j*n_grid+i] = (i+0.5)*g_size-length;
      posy[j*n_grid+i] = (j+0.5)*g_size-length;
    }
  }

  for(i=0;i<cones[0];i++){
    srand(time(NULL));
    c_z[0] = sin(cones[1+4*i])*cos(cones[2+4*i]);
    c_z[1] = sin(cones[1+4*i])*sin(cones[2+4*i]);
    c_z[2] = cos(cones[1+4*i]);
    vec_ort_wc(c_x,d_x,c_z);
    vec_cross_wc(c_y,c_x,c_z);
    vec_rotate(c_y,PI*0.5-cones[4+4*i],c_x);
    for(j=0;j<10000;j++){
      phi = (double) rand()/RAND_MAX;
      phi *= 2*PI;
      vec_rotate_wc(dir,c_z,phi,c_x);
      ratio = height/dir[2];
      dir[0] *= ratio;
      dir[1] *= ratio;
      if(dir[2] > 0 && dir[0] < length && dir[0] > -length && dir[1] < length && dir[1] > -length){
	for(k=0;k<n_grid*n_grid;k++){
	  if(dir[0] <= posx[k]+g_size_over_2 && dir[0] > posx[k]-g_size_over_2 &&
	     dir[1] <= posy[k]+g_size_over_2 && dir[1] > posy[k]-g_size_over_2){
	    image[k] += cones[3+4*i]/10000; break;
	  }
	}
      }
    }
  }

  FILE* f = fopen("image.txt","w+");
  for(i=0;i<n_grid;i++){
    for(j=0;j<n_grid;j++){
      fprintf(f,"%f ",image[j*n_grid+i]);
    }
    fprintf(f,"\n");
  }
  fclose(f);

  f = fopen("image.tpt","w+");
  for(i=0;i<n_grid;i++){
    for(j=0;j<n_grid;j++){
      fprintf(f,"%f %f %f \n",(float)i,(float)j,image[j*n_grid+i]);
    }
  }
  fclose(f);

  return 0;
}

/*****************************************************************************/
