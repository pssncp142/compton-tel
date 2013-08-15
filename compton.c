/*****************************************************************************\
 * Yigit Dallilar 10.07.2013                                                 *
 * DTU-Space                                                                 *
 * Compton telescope looks for interaction order                             *
\*****************************************************************************
 * Variables :                                                               *
 * path[0]                         --> length of a path                      *
 * path[1]                         --> number of paths                       *
 * path[2+path[0]*i+0:(path[0]-1)] --> interaction list                      *
 * i                               --> path index                            *
 *                                                                           *
 * match[0]                        --> number of matches                     *
 * match[1+5*i+3]                  --> angle calculated from positions       *
 * match[1+5*i+4]                  --> angle calculated from compton         *
 * if first match for the event                                              *
 * match[1+5*i+0:2]                --> direction from absorption to compton  *
 * for other matches                                                         *
 * match[1+5*i]                    --> path index                            *
 * match[1+5*i+1]                  --> last point for the current path       *
 * match[1+5*i+2]                  --> the point investigated                *
 *****************************************************************************/

#include "stdio.h"
#include "math.h"
#include "compton.h"
#include "3d_cart_vec.h"
#include "event.h"

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

#define PI 3.14159f

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

static int i,j,k,l,m,n;
static double rest_e = 510.998;

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//returns compton angle for the energies
double compt_angle(double post, double pre){
  
  return acos(1+rest_e*(post-pre)/(pre*post));
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

double compt_en_err_en(double ang, double post, double err_en){

  return err_en/(pow(1/compt_com_en(ang,post),2)*post);
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

double compt_en_err_pos1(double ang, double post, double dist, double err_pos){

  return err_pos*sin(ang)/(pow(1/compt_com_en(ang,post),2)*post*dist);
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

double compt_en_err_pos2(double ang, double post, double dist, double err_pos){

  double err_dist=err_pos/dist;

  return (sin(ang)*err_dist+0.75*pow(sin(ang)*err_dist,2)+0.5*cos(ang)*err_dist*err_dist)/
    (pow(1/compt_com_en(ang,post),2)*post);
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

double compt_com_en(double ang, double post){

  return 1/(1/post-(1-cos(ang))/rest_e);
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//
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
  }

  if(l_p!=0){
    if(verb) printf("Path  It Post    CompDep(keV)    DetDep(keV)    AbsDiff    DetDep0.2    PosAng    ComptAng    ErrAng\n");
    for(i=0;i<n_p;i++){
      pre_e=0;
      for(j=0;j<l_p-1;j++) pre_e += event[1+4*path[2+l_p*i+j]];
      for(j=0;j<n_ev;j++){
	fnd = 1;
	for(k=0;k<l_p;k++){
	  if(j==path[2+l_p*i+k]) {fnd = 0; break;}
	}
	if(fnd){
	  for(k=0;k<3;k++){
	    vec1[k] = event[2+4*path[2+l_p*i+l_p-1]+k] - event[2+4*j+k];
	    vec2[k] = event[2+4*path[2+l_p*i+l_p-2]+k] - event[2+4*path[2+l_p*i+l_p-1]+k];
	  }
	  match[ndx++] = i;
	  match[ndx++] = path[2+l_p*i+l_p-1];
	  match[ndx++] = j;
	  match[ndx++] = compt_com_en(acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2))),pre_e)-pre_e;
	  match[ndx++] = event[1+4*path[2+l_p*i+l_p-1]];
	  match[ndx++] = acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2)));
	  match[ndx++] = compt_angle(pre_e,pre_e+event[1+4*path[2+l_p*i+l_p-1]]);
	  match[0]++;
	  if(verb) printf("%4d %3d %4d %15.3f %14.3f %10.3f %12.3f %9.3f %11.3f %9.3f\n",
			  (int)match[ndx-7],(int)match[ndx-6],(int)match[ndx-5],match[ndx-4],match[ndx-3],
			  fabs(match[ndx-4]-match[ndx-3]),match[ndx-3]*0.2,match[ndx-2],match[ndx-1],fabs(match[ndx-2]-match[ndx-1]));
	}
      }
    }
  } else {
    if(verb) printf("Pre  It Post    CompDep(keV)    DetDep(keV)    AbsDiff    DetDep0.2    PosAng    ComptAng    ErrAng\n");
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
	    match[ndx++] = compt_com_en(acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2))),event[1+4*k])-event[1+4*k];
	    match[ndx++] = event[1+4*j];	  
	    match[ndx++] = acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2)));
	    match[ndx++] = compt_angle(event[1+4*k],event[1+4*k]+event[1+4*j]);	  
	    if(verb) printf("%3d %3d %4d %15.3f %14.3f %10.3f %12.3f %9.3f %11.3f %9.3f\n",
			    (int)match[ndx-7],(int)match[ndx-6],(int)match[ndx-5],match[ndx-4],match[ndx-3],
			    fabs(match[ndx-4]-match[ndx-3]),match[ndx-3]*0.2,match[ndx-2],match[ndx-1],fabs(match[ndx-2]-match[ndx-1]));
	    
	  }
	}
      }
    }
    
  }

  return 0;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int compt_add_path(int path[], double match[], double event[], int verb){
  
  int s = 0;
  int cnt = 0;
  int tpath[10000] = {0};
  int ndx = 0;
  int ndx_p,fnd;
  int n_p=0;

  if(path[0] == 0){
    for(i=0;i<match[0];i++){
      if(fabs(match[1+7*i+3]-match[1+7*i+4])<0.2*match[1+7*i+4]){
	path[0]=2;
	s = 1;
	path[1]++;
	path[2+3*cnt] = (int)match[1+7*i];
	path[2+3*cnt+1] = (int)match[1+7*i+1];
	path[2+3*cnt+2] = (int)match[1+7*i+2];
	cnt++;
      }
    }
  } else if (match[0]==0){
    
  } else {
    for(i=0;i<match[0];i++){
      if(fabs(match[1+7*i+3]-match[1+7*i+4])<0.2*match[1+7*i+4]){
	for(j=0;j<path[1];j++){
	  if(j == (int)match[1+7*i]){
	    fnd=1;
	    if(fnd){
	      s=1; n_p++;
	      for(k=0;k<path[0];k++){
		tpath[ndx+k] = path[2+path[0]*j+k];
	      }
	      tpath[ndx+path[0]] = (int)match[1+7*i+2];
	      cnt++;
	      ndx += (path[0]+1);
	    }
	  }
	}
      }
    }
    path[1] = n_p;
    for(j=0;j<5000;j++) {path[2+j] = tpath[j];}
  }

  if(verb){
    printf("Real interaction order :\n");
    printf("* ");
    if(event[1+(int)event[0]*5]==0) printf(" Out -->");
    for(j=0;j<event[0];j++){
      printf("%2d -->",(int)event[1+(int)event[0]*4+j]);
    }
    printf(" Sky\n");
  }   

  if(s){
    path[0]++;
    if(verb){
      printf("Number of paths : %2d\n",path[1]);
      for(i=0;i<path[1];i++){
	printf("* ");
	for(j=0;j<path[0];j++){
	  printf("%2d -->",path[2+path[0]*i+j]);
	}
	printf(" ?\n");
      }
    }
    return 0;
  } else {
    if(verb){
      if(path[1]==0){
	printf("\n\nUncompleted event...\n");
      } else {
	printf("\n\nResulted paths : %2d\n",path[1]);
	for(i=0;i<path[1];i++){
	  printf("* ");
	  for(j=0;j<path[0];j++){
	    printf("%2d -->",path[2+path[0]*i+j]);
	  }
	  printf(" Sky\n");
	}
      }
    }
    return 1;
  }
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//gets the full path

int compt_get_path(int path[], double event[], int verb){

  double match[10000];

  path[0] = 0; path[1]=0;
  for(j=0;j<100;j++){
    compt_match3(match,path,event,verb);
    if(compt_add_path(path,match,event,verb)==1) break;
  }

  if(event[1+5*(int)event[0]]==1){
    if(path[1]==1){
      return 0;
    } else if(path[1]==0){
      return 1;
    } else {
      return 4;
    }
  } else {
    if(path[1]==0){
      return 2;
    } else {
      return 3;
    }
  }
}

/*****************************************************************************/
