/*****************************************************************************\
 * Yigit Dallilar 10.07.2013                                                 *
 * DTU-Space                                                                 *
 * Compton telescope returns data with specific parameters defined.          *
\*****************************************************************************/

#include "stdio.h"
#include "exclude.h"

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int f_ndx,l_ndx,e_ndx;
int i,j,k;
int tot_ev;
int fnd,cnt;
const int v_l = 7;
const int p_l = 13;

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int energy_int(double data[], int n_proc[], int *n_ev, double en[]){

  tot_ev = *n_ev;
  cnt = 0;
  f_ndx=0;
  e_ndx=0;
  for(i=0;i<tot_ev;i++) e_ndx += v_l + n_proc[i]*p_l; 

  for(i=0;i<tot_ev;i++){
    fnd = 0;
    l_ndx = f_ndx + v_l + n_proc[i-cnt]*p_l;
    if(data[f_ndx+6]<en[0] || data[f_ndx+6]>en[1]){
      fnd = 1;
    }
    if(fnd){
      for(j=0;j<e_ndx-l_ndx;j++) data[f_ndx+j] = data[l_ndx+j];
      e_ndx -= (v_l+n_proc[i-cnt]*p_l);
      for(j=0;j<tot_ev-i;j++) n_proc[i-cnt+j] = n_proc[i-cnt+j+1];
      cnt++;
    } else {
      f_ndx = l_ndx;
    }
  }

  *n_ev = tot_ev - cnt;

  return 0;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int nofproc(double data[], int n_proc[], int *n_ev, int n[]){

  tot_ev = *n_ev;
  cnt = 0;
  f_ndx=0;
  e_ndx=0;
  for(i=0;i<tot_ev;i++) e_ndx += v_l + n_proc[i]*p_l; 

  for(i=0;i<tot_ev;i++){
    fnd = 0;
    l_ndx = f_ndx + v_l + n_proc[i-cnt]*p_l;
    if(n_proc[i-cnt]<n[0] || n_proc[i-cnt]>n[1]){
      fnd = 1;
    }
    if(fnd){
      for(j=0;j<e_ndx-l_ndx;j++) data[f_ndx+j] = data[l_ndx+j];
      e_ndx -= (v_l+n_proc[i-cnt]*p_l);
      for(j=0;j<tot_ev-i;j++) n_proc[i-cnt+j] = n_proc[i-cnt+j+1];
      cnt++;
    } else {
      f_ndx = l_ndx;
    }
  }

  *n_ev = tot_ev - cnt;

  return 0;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int clean_0(double data[], int n_proc[], int *n_ev){

  tot_ev = *n_ev;
  cnt = 0;
  f_ndx=0;
  e_ndx=0;
  for(i=0;i<tot_ev;i++) e_ndx += v_l + n_proc[i]*p_l; 

  for(i=0;i<tot_ev;i++){
    fnd = 0;
    l_ndx = f_ndx + v_l + n_proc[i-cnt]*p_l;
    if(n_proc[i-cnt]==0){
      fnd = 1;
    }
    if(fnd){
      for(j=0;j<e_ndx-l_ndx;j++) data[f_ndx+j] = data[l_ndx+j];
      e_ndx -= (v_l+n_proc[i-cnt]*p_l);
      for(j=0;j<tot_ev-i;j++) n_proc[i-cnt+j] = n_proc[i-cnt+j+1];
      cnt++;
    } else {
      f_ndx = l_ndx;
    }
  }

  *n_ev = tot_ev - cnt;
  return 0;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int proc_exc(double data[], int n_proc[], int *n_ev, char proc[]){
  
  char c_proc[3] = {'c','p','r'};
  double check = 3;

  tot_ev = *n_ev;
  cnt=0;
  f_ndx = 0;  
  e_ndx = 0;
  for(i=0;i<tot_ev;i++) e_ndx += v_l + n_proc[i]*p_l; 

  if (proc[0] == c_proc[0]){
    check = 0;
  } else if (proc[0] == c_proc[1]){
    check = 1;
  } else if (proc[0] == c_proc[2]){
    check = 2;
  } else {
    printf("!!! Invalid process...\n");
    goto exit;
  }

  for(i=0;i<tot_ev;i++){
    fnd = 0;
    l_ndx = f_ndx + v_l + n_proc[i-cnt]*p_l;
    for(j=0;j<n_proc[i-cnt];j++){
      if(data[f_ndx+v_l+p_l*j] == check){
	fnd = 1;
	break;
      }
    }
    if(fnd){
      for(j=0;j<e_ndx-l_ndx;j++) data[f_ndx+j] = data[l_ndx+j];
      e_ndx -= (v_l+n_proc[i-cnt]*p_l);
      for(j=0;j<tot_ev-i;j++) n_proc[i-cnt+j] = n_proc[i-cnt+j+1];
      cnt++;
    } else {
      f_ndx = l_ndx;
    }
  }

  *n_ev = tot_ev - cnt;

 exit :
  if (check==3)
    return 1;
  else
    return 0;

}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int all_detected(double data[], int n_proc[], int *n_ev){
  

  tot_ev = *n_ev;
  cnt=0;
  f_ndx = 0;  
  e_ndx = 0;
  double ener;

  for(i=0;i<tot_ev;i++) e_ndx += v_l + n_proc[i]*p_l; 
  
  for(i=0;i<tot_ev;i++){
    fnd = 0;
    ener = 0;
    l_ndx = f_ndx + v_l + n_proc[i-cnt]*p_l;
    for(j=0;j<n_proc[i-cnt];j++){
      ener += data[f_ndx+v_l+p_l*j+12];
    }
    if(data[f_ndx+6] != ener) fnd=1;
    if(fnd){
      for(j=0;j<e_ndx-l_ndx;j++) data[f_ndx+j] = data[l_ndx+j];
      e_ndx -= (v_l+n_proc[i-cnt]*p_l);
      for(j=0;j<tot_ev-i;j++) n_proc[i-cnt+j] = n_proc[i-cnt+j+1];
      cnt++;
    } else {
      f_ndx = l_ndx;
    }
  }
  
  *n_ev = tot_ev - cnt;

  return 0;

}

/*****************************************************************************/