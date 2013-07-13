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
int v_l;
int p_l;

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//exclude event from data within energy range considering en array
int energy_int(double data[], int n_proc[], int *n_ev, double en[]){

  v_l =  7;
  p_l = 13;

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

  printf("Data option enabled -> in the range [%5.2f,%f5.2] keV...\n",en[0],en[1]);

  return 0;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//exclude events from data within the range of n (# of interaction)
int nofproc(double data[], int n_proc[], int *n_ev, int n[]){

  v_l =  7;
  p_l = 13;

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

  printf("Data option enabled -> # of process [%2d,%2d] ...\n",n[0],n[1]);

  return 0;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//exclude all events with no interactions
int clean_0(double data[], int n_proc[], int *n_ev){

  v_l =  7;
  p_l = 13;

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

  printf("Data option enabled -> Neglect 0 interaction...\n");

  return 0;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//exclude events with specific interaction
int proc_exc(double data[], int n_proc[], int *n_ev, char proc[]){

  v_l =  7;
  p_l = 13;
  
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

  printf("Data option enabled -> exclude process %s ...\n",proc);

 exit :
  if (check==3)
    return 1;
  else
    return 0;

}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//to get only complete events
int all_detected(double data[], int n_proc[], int *n_ev){

  v_l =  7;
  p_l = 13;
  
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

  printf("Data option enabled -> Only events detected completely...\n");


  return 0;

}

/*****************************************************************************/
