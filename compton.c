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
#include "stdlib.h"

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
  int *tpath = (int*)malloc(1000000*sizeof(int));
  int ndx = 0;
  int ndx_p,fnd;
  int n_p=0;

  if(path[0] == 0){
    for(i=0;i<match[0];i++){
      if(fabs(match[1+7*i+3]-match[1+7*i+4])<0.5*match[1+7*i+4]){
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
    for(j=0;j<1000000;j++) {path[2+j] = tpath[j];}
  }

  free(tpath);

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

int compt_analyse(double *event, int *path){

  double cone[15];

  int succ=0;
  int i,j,k;
  int type = event_type(event);
  cone[0] = (double) type;
  //printf("%d\n",type);

  if(type==2){

    double en;
    int check=0;
    int other;
    for(i=0;i<2;i++){
      other=!i;
      en = event[1+4*other]*rest_e/(rest_e+2*event[1+4*other]);
      if(event[1+4*i]>en) check += i+1;
    }
    check--;
    cone[1] = (double)check;
    cone[2] = event[2+4*check];
    cone[3] = event[3+4*check];
    cone[4] = event[4+4*check];
    cone[5] = event[2+4*(!check)];
    cone[6] = event[3+4*(!check)];
    cone[7] = event[4+4*(!check)];
    cone[8] = event[1+4*check]+event[1+4*(!check)];
    succ=0;

  } else if(type==3){

    double post,pre;
    double vec1[3],vec2[3];
    double err[6];
    int comb[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}}; 
    double en1,en2;
    double ang;
    for(i=0;i<6;i++){
      vec1[0] = event[2+4*comb[i][1]] - event[2+4*comb[i][2]];  
      vec1[1] = event[3+4*comb[i][1]] - event[3+4*comb[i][2]];  
      vec1[2] = event[4+4*comb[i][1]] - event[4+4*comb[i][2]];  
      vec2[0] = event[2+4*comb[i][0]] - event[2+4*comb[i][1]];  
      vec2[1] = event[3+4*comb[i][0]] - event[3+4*comb[i][1]];  
      vec2[2] = event[4+4*comb[i][0]] - event[4+4*comb[i][1]];
      en1 = event[1+4*comb[i][0]];
      en2 = event[1+4*comb[i][0]]+event[1+4*comb[i][1]];
      ang = acos(1+rest_e*(en1-en2)/(en1*en2));
      err[i] = fabs(ang-acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2))));
    }
    double min=1e6;
    int ndx;
    for(i=0;i<6;i++){
      if(min>err[i]) {min = err[i]; ndx=i;}
    }
    if(ndx>=0 && ndx<6){
      cone[1] = event[2+4*comb[ndx][2]];
      cone[2] = event[3+4*comb[ndx][2]];
      cone[3] = event[4+4*comb[ndx][2]];
      cone[4] = event[2+4*comb[ndx][1]];
      cone[5] = event[3+4*comb[ndx][1]];
      cone[6] = event[4+4*comb[ndx][1]];
      cone[7] = event[2+4*comb[ndx][0]];
      cone[8] = event[3+4*comb[ndx][0]];
      cone[9] = event[4+4*comb[ndx][0]];
      cone[10] = event[1+4*comb[ndx][0]]+event[1+4*comb[ndx][1]]+event[1+4*comb[ndx][2]];
      cone[11] = err[ndx];
      succ=0;
    }
  } /*else if (type==9){

    printf("a\n");
    int i,j,k;
    double x=25,y=25,z=-50;
    double vec1[3],vec2[3];
    int ins;
    int out[10];
    int cnt1=0;
    int cnt2=0;
    int st;
    double ev[200];
    double en;

    for(i=0;i<event[0];i++){
      if(event[2+4*i]<x && event[2+4*i]>-x &&
	 event[3+4*i]<y && event[3+4*i]>-y &&
	 event[4+4*i]>z){
	ins=i;
	cnt1++;
      } else {
	out[cnt2]=i;
	cnt2++;
      }
    }
    
    ev[0]=(double)cnt2;
    printf("cnt : %d\n",cnt2);
    for(i=0;i<cnt2;i++){
      ev[1+4*i] = event[1+4*out[i]];
      ev[2+4*i] = event[2+4*out[i]];
      ev[3+4*i] = event[3+4*out[i]];
      ev[4+4*i] = event[4+4*out[i]];
      printf("A %f\n",ev[1+4*i]);
    }
    st=compt_get_path(path,ev,0);
    printf("cnt : %d %f %d\n",cnt2,ev[0],path[1]);

    en=0;
    for(j=0;j<path[0];j++){
      en += ev[1+4*path[2+j]];
    } 
    printf("%f\n",en);


    double en1,en2;
    double ang;
    double err[path[1]];
   
    if(path[1]!=0){
      //printf("%d %d\n",path[0],path[1]);
      for(i=0;i<path[1];i++){
	en=0;
	for(j=0;j<path[0]-1;j++){
	  en += ev[1+4*path[2+j+path[0]*i]];
	}
	printf("%f\n",en);
	vec1[0] = ev[2+4*path[1+path[0]*(i+1)]] - event[2+4*ins];  
	vec1[1] = ev[3+4*path[1+path[0]*(i+1)]] - event[3+4*ins];  
	vec1[2] = ev[4+4*path[1+path[0]*(i+1)]] - event[4+4*ins];  
	vec2[0] = ev[2+4*path[path[0]*(i+1)]] - ev[2+4*path[1+path[0]*(i+1)]];  
	vec2[1] = ev[3+4*path[path[0]*(i+1)]] - ev[3+4*path[1+path[0]*(i+1)]];  
	vec2[2] = ev[4+4*path[path[0]*(i+1)]] - ev[4+4*path[1+path[0]*(i+1)]];
	en1 = en;
	en2 = en+ev[1+4*path[1+path[0]*(i+1)]];
	ang = acos(1+rest_e*(en1-en2)/(en1*en2));
	err[i] = fabs(ang-acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2))));      
      }
            
      double min=1e6;
      int ndx=80;
      for(i=0;i<path[1];i++){
	if(min>err[i]) { min = err[i]; ndx=i;}
	printf("%f %d %d\n",err[i],i,path[1]);
      }
      printf("%d\n",ndx);
      if(ndx>=0 && ndx<path[1]){
	cone[1] = event[2+4*ins];
	cone[2] = event[3+4*ins];
	cone[3] = event[4+4*ins];
	cone[4] = ev[2+4*path[1+path[0]*(ndx+1)]];
	cone[5] = ev[3+4*path[1+path[0]*(ndx+1)]];
	cone[6] = ev[4+4*path[1+path[0]*(ndx+1)]];
	cone[7] = ev[2+4*path[1+path[0]*(ndx+1)]];
	cone[8] = ev[2+4*path[1+path[0]*(ndx+1)]];
	cone[9] = ev[2+4*path[1+path[0]*(ndx+1)]];
	en=0;
	for(j=0;j<path[0];j++){
	  en += ev[1+4*path[2+j+path[1]*ndx]];
	}
	cone[10] = en+event[1+4*ins];
	cone[11] = err[ndx];
	succ=1;
      } 
    }
    
    }*/ else if(type==5){

    double x=25,y=25,z=-50;
    int ins;
    for(i=0;i<2;i++){
      if(event[2+4*i]<x && event[2+4*i]>-x &&
	 event[3+4*i]<y && event[3+4*i]>-y &&
	 event[4+4*i]>z){
	ins = i;
      }
    }

    cone[1] = (double)ins;
    cone[2] = event[2+4*ins];
    cone[3] = event[3+4*ins];
    cone[4] = event[4+4*ins];
    cone[5] = event[2+4*(!ins)];
    cone[6] = event[3+4*(!ins)];
    cone[7] = event[4+4*(!ins)];
    cone[8] = event[1+4*ins]+event[1+4*(!ins)];
    succ=1;
    
  } else if(type==6){
    
    double x=25,y=25,z=-50;
    double vec1[3],vec2[3];
    int ins[2];
    int out;
    int cnt=0;
    for(i=0;i<3;i++){
      if(event[2+4*i]<x && event[2+4*i]>-x &&
	 event[3+4*i]<y && event[3+4*i]>-y &&
	 event[4+4*i]>z){
	ins[cnt]=i;
	cnt++;
      } else {
	out = i;
      }
    }
    
    double en1,en2;
    double ang;
    double err[2];
    for(i=0;i<2;i++){
      vec1[0] = event[2+4*ins[!i]] - event[2+4*ins[i]];  
      vec1[1] = event[3+4*ins[!i]] - event[3+4*ins[i]];  
      vec1[2] = event[4+4*ins[!i]] - event[4+4*ins[i]];  
      vec2[0] = event[2+4*out] - event[2+4*ins[!i]];  
      vec2[1] = event[3+4*out] - event[3+4*ins[!i]];  
      vec2[2] = event[4+4*out] - event[4+4*ins[!i]];
      en1 = event[1+4*out];
      en2 = event[1+4*out]+event[1+4*ins[!i]];
      ang = acos(1+rest_e*(en1-en2)/(en1*en2));
      err[i] = fabs(ang-acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2))));      
    }

    double min=1e6;
    int ndx=5;
    for(i=0;i<2;i++){
      if(min>err[i]) {min = err[i]; ndx=i;}
    }
    if(ndx==0 || ndx==1){
      cone[1] = event[2+4*ins[ndx]];
      cone[2] = event[3+4*ins[ndx]];
      cone[3] = event[4+4*ins[ndx]];
      cone[4] = event[2+4*ins[!ndx]];
      cone[5] = event[3+4*ins[!ndx]];
      cone[6] = event[4+4*ins[!ndx]];
      cone[7] = event[2+4*out];
      cone[8] = event[3+4*out];
      cone[9] = event[4+4*out];
      cone[10] = event[1+4*ins[ndx]]+event[1+4*ins[!ndx]]+event[1+4*out];
      cone[11] = err[ndx];
      succ=1;
    } 
    
  } else if(type==8){

    double x=25,y=25,z=-50;
    double vec1[3],vec2[3];
    int ins;
    int out[2];
    int cnt=0;
    for(i=0;i<3;i++){
      if(event[2+4*i]<x && event[2+4*i]>-x &&
	 event[3+4*i]<y && event[3+4*i]>-y &&
	 event[4+4*i]>z){
	ins=i;
      } else {
	out[cnt] = i;
	cnt++;
      }
    }

    double en1,en2;
    double ang;
    double err[2];
    for(i=0;i<2;i++){
      vec1[0] = event[2+4*out[i]] - event[2+4*ins];  
      vec1[1] = event[3+4*out[i]] - event[3+4*ins];  
      vec1[2] = event[4+4*out[i]] - event[4+4*ins];  
      vec2[0] = event[2+4*out[!i]] - event[2+4*out[i]];  
      vec2[1] = event[3+4*out[!i]] - event[3+4*out[i]];  
      vec2[2] = event[4+4*out[!i]] - event[4+4*out[i]];
      en1 = event[1+4*out[!i]];
      en2 = event[1+4*out[!i]]+event[1+4*out[i]];
      ang = acos(1+rest_e*(en1-en2)/(en1*en2));
      err[i] = fabs(ang-acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2))));      
    }

    double min=1e6;
    int ndx;
    for(i=0;i<2;i++){
      if(min>err[i]) {min = err[i]; ndx=i;}
    }
    if(ndx==0 || ndx==1){
      cone[1] = event[2+4*ins];
      cone[2] = event[3+4*ins];
      cone[3] = event[4+4*ins];
      cone[4] = event[2+4*out[ndx]];
      cone[5] = event[3+4*out[ndx]];
      cone[6] = event[4+4*out[ndx]];
      cone[7] = event[2+4*out[!ndx]];
      cone[8] = event[3+4*out[!ndx]];
      cone[9] = event[4+4*out[!ndx]];
      cone[10] = event[1+4*ins]+event[1+4*out[ndx]]+event[1+4*out[!ndx]];
      cone[11] = err[ndx];
      succ=1;
    }

  } else if (type==9){
    
    double x=25,y=25,z=-50;
    double vec1[3],vec2[3];
    int ins;
    int out[20];
    int cnt2=0;

    for(i=0;i<event[0];i++){
      if(event[2+4*i]<x && event[2+4*i]>-x &&
	 event[3+4*i]<y && event[3+4*i]>-y &&
	 event[4+4*i]>z){
	ins=i;
      } else {
	out[cnt2]=i;
	cnt2++;
      }
    }

    double en1,en2;
    double en=0;
    double ang;
    double err[cnt2*cnt2];
    for(i=0;i<cnt2;i++){
      for(j=0;j<cnt2;j++){
	if(j!=k){
	  en=0;
	  for(k=0;k<cnt2;k++){
	    if(k!=j)
	      en += event[1+4*out[k]];
	  }
	  vec1[0] = event[2+4*out[j]] - event[2+4*ins];  
	  vec1[1] = event[3+4*out[j]] - event[3+4*ins];  
	  vec1[2] = event[4+4*out[j]] - event[4+4*ins];  
	  vec2[0] = event[2+4*out[i]] - event[2+4*out[j]];  
	  vec2[1] = event[3+4*out[i]] - event[3+4*out[j]];  
	  vec2[2] = event[4+4*out[i]] - event[4+4*out[j]];
	  en1 = en;
	  en2 = en + event[1+4*out[j]];
	  ang = acos(1+rest_e*(en1-en2)/(en1*en2));
	  err[i*cnt2+j] = fabs(ang-acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2))));      
	}
      }
    }

    double min=1e6;
    int ndx[2]={50,50};
    for(i=0;i<cnt2;i++){
      for(j=0;j<cnt2;j++){
	if(min>err[i*cnt2+j] && i!=j) {min = err[i*cnt2+j]; ndx[0]=i; ndx[1]=j;}
      }
    }

    if(ndx[0]>=0 && ndx[1]>=0 && ndx[0]<cnt2 && ndx[1]<cnt2){
      cone[1] = event[2+4*ins];
      cone[2] = event[3+4*ins];
      cone[3] = event[4+4*ins];
      cone[4] = event[2+4*out[ndx[1]]];
      cone[5] = event[3+4*out[ndx[1]]];
      cone[6] = event[4+4*out[ndx[1]]];
      cone[7] = event[2+4*out[ndx[0]]];
      cone[8] = event[3+4*out[ndx[0]]];
      cone[9] = event[4+4*out[ndx[0]]];
      en=0;
      for(k=0;k<cnt2;k++){
	  en += event[1+4*out[k]];
      }
      cone[10] = en + event[1+4*ins];
      cone[11] = err[ndx[0]*cnt2+ndx[1]];
      succ=1;
    }


  } else if (type == 10){
    
    double x=25,y=25,z=-50;
    double vec1[3],vec2[3];
    int ins[2];
    int out[2];
    int cnt1=0;
    int cnt2=0;

    for(i=0;i<4;i++){
      if(event[2+4*i]<x && event[2+4*i]>-x &&
	 event[3+4*i]<y && event[3+4*i]>-y &&
	 event[4+4*i]>z){
	ins[cnt1]=i;
	cnt1++;
      } else {
	out[cnt2]=i;
	cnt2++;
      }
    }

    double en1,en2;
    double ang;
    double err[4];
    for(i=0;i<2;i++){
      for(j=0;j<2;j++){
	vec1[0] = event[2+4*ins[!i]] - event[2+4*ins[i]];  
	vec1[1] = event[3+4*ins[!i]] - event[3+4*ins[i]];  
	vec1[2] = event[4+4*ins[!i]] - event[4+4*ins[i]];  
	vec2[0] = event[2+4*out[j]] - event[2+4*ins[!i]];  
	vec2[1] = event[3+4*out[j]] - event[3+4*ins[!i]];  
	vec2[2] = event[4+4*out[j]] - event[4+4*ins[!i]];
	en1 = event[1+4*out[!j]]+event[1+4*out[j]];
	en2 = event[1+4*out[!j]]+event[1+4*out[j]]+event[1+4*ins[!i]];
	ang = acos(1+rest_e*(en1-en2)/(en1*en2));
	err[i*2+j] = fabs(ang-acos(vec_dotp(vec1,vec2)/(vec_norm(vec1)*vec_norm(vec2))));      
      }
    }
    
    double min=1e6;
    int ndx[2];
    for(i=0;i<2;i++){
      for(j=0;j<2;j++){
	if(min>err[i*2+j]) {min = err[i*2+j]; ndx[0]=i; ndx[1]=j;}
      }
    }

    if(ndx[0]==0 || ndx[0]==1 || ndx[1]==0 || ndx[1]==1){
      cone[1] = event[2+4*ins[ndx[0]]];
      cone[2] = event[3+4*ins[ndx[0]]];
      cone[3] = event[4+4*ins[ndx[0]]];
      cone[4] = event[2+4*ins[!ndx[0]]];
      cone[5] = event[3+4*ins[!ndx[0]]];
      cone[6] = event[4+4*ins[!ndx[0]]];
      cone[7] = event[2+4*out[ndx[1]]];
      cone[8] = event[3+4*out[ndx[1]]];
      cone[9] = event[4+4*out[ndx[1]]];
      cone[10] = event[1+4*ins[ndx[0]]]+event[1+4*ins[!ndx[0]]]+event[1+4*out[ndx[1]]]+event[1+4*out[!ndx[1]]];
      cone[11] = err[ndx[0]*2+ndx[1]];
      succ=1;
    }

  }

  if(succ)
    printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
	   cone[0],cone[1],cone[2],cone[3],cone[4],cone[5],cone[6],cone[7],cone[8],cone[9],cone[10],cone[11]);

  return 0;
}

/*****************************************************************************/
