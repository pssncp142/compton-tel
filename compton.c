/*****************************************************************************\
 * Yigit Dallilar 10.07.2013                                                 *
 * DTU-Space                                                                 *
 * Compton telescope looks for interaction order                             *
\*****************************************************************************/

#include "stdio.h"
#include "math.h"
#include "compton.h"
#include "3d_cart_vec.h"

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

#define PI 3.14159f

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

static int i,j,k,l,m,n;
static double rest_e = 510.998;

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

double compt_angle(double post, double pre){
  
  return PI-acos(rest_e*(pre-post)/(pre*post)-1);
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int compt_match3(double match[], int path[], double event[], int verb){

  int l_p,n_p,n_ev;
  int fnd,ndx;
  double pre_e;
  double vec1[5];
  double vec2[5];
  l_p=path[0];
  n_p=path[1];
  n_ev=event[0];
  match[0]=0;
  ndx=1;

  if (verb){
    printf("\nMatching results :\n");
    printf("Pre  It Post PosAng ComptAng  Error\n");
  }

  if(l_p!=0){
    for(i=0;i<n_p;i++){
      pre_e=0;
      for(j=0;j<l_p-1;j++) pre_e += event[1+4*path[2+l_p*i+j]];
      for(j=0;j<n_ev;j++){
	fnd = 0;
	for(k=0;k<l_p;k++){
	  if(j==path[2+l_p*i+k]) {fnd = 1; break;}
	}
	if(!fnd){
	  for(k=0;k<3;k++){
	    vec1[k] = event[2+4*path[2+l_p*i+l_p-1]+k] - event[2+4*j+k];
	    vec2[k] = event[2+4*path[2+l_p*i+l_p-2]+k] - event[2+4*path[2+l_p*i+l_p-1]+k];
	  }
	  match[ndx++] = path[2+l_p*i+l_p-2];
	  match[ndx++] = path[2+l_p*i+l_p-1];
	  match[ndx++] = j;
	  match[ndx++] = acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2)));
	  match[ndx++] = compt_angle(pre_e,pre_e+event[1+4*path[2+l_p*i+l_p-1]]);
	  if(i==0) match[0]++;
	  if(verb) printf("%3d %3d %4d %6.3f %8.3f %6.3f\n",
			  (int)match[ndx-5],(int)match[ndx-4],(int)match[ndx-3],match[ndx-2],match[ndx-1],match[ndx-2]-match[ndx-1]);
	}
      }
    }
  } else {
    for(j=0;j<n_ev;j++){
      for(k=0;k<n_ev;k++){
	for(l=0;l<n_ev;l++){
	  if(j!=k && k!=l && j!=l){
	    match[0]++;
	    for(m=0;m<3;m++){
	      vec1[m] = event[2+4*j+m] - event[2+4*l+m];
	      vec2[m] = event[2+4*k+m] - event[2+4*j+m];
	    }
	    match[ndx++] = k;
	    match[ndx++] = j;
	    match[ndx++] = l;
	    match[ndx++] = acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2)));
	    match[ndx++] = compt_angle(event[1+4*k],event[1+4*k]+event[1+4*j]);	    
	    if(verb) printf("%4d %3d %3d %6.3f %8.3f %6.3f\n",
			    (int)match[ndx-5],(int)match[ndx-4],(int)match[ndx-3],match[ndx-2],match[ndx-1],match[ndx-2]-match[ndx-1]);
	  }
	}
      }
    }
    
  }

  return 0;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int compt_add_path(int path[], double match[], int verb){
  
  int s = 0;
  int cnt = 0;

  if(path[0] == 0){
    for(i=0;i<match[0];i++){
      if(fabs(match[1+5*i+3]-match[1+5*i+4])<0.5){
	path[0]=2;
	s = 1;
	path[1]++;
	path[2+3*cnt] = (int)match[1+5*i];
	path[2+3*cnt+1] = (int)match[1+5*i+1];
	path[2+3*cnt+2] = (int)match[1+5*i+2];
	cnt++;
      }
    }
  } else {
    
  }
 
  if(s){
    path[0]++;
    if(verb){
      printf("Number of paths : %2d\n",path[1]);
      for(i=0;i<path[1];i++){
	printf(" * ");
	for(j=0;j<path[0];j++){
	  printf("%3d -->",path[2+path[0]*i+j]);
	}
	printf("?\n");
      }
    }
    return 0;
  } else {
    if(verb) printf("No track added...\n");
    return 1;
  }
}

/*****************************************************************************/
