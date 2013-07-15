/*****************************************************************************\
 * Yigit Dallilar 12.07.2013                                                 *
 * DTU-Space                                                                 *
 * Compton telescope library for events                                      *
\*****************************************************************************
 * Variables :                                                               *
 * data   --> direct input from geant binary                                 *
 * event  --> event array                                                    *
 * n_proc --> # of process related with data                                 *
 * n_ev   --> from # of events                                               *
 * ndx    --> event index                                                    *
 * verb   --> verbosity                                                      *
 *                                                                           *
 * event[0]          --> number of interaction                               *
 * event[1+4*j]      --> interaction energy                                  *
 * event[2+4*j+0:2]  --> interaction positions                               *
 * end               --> last index for the interactions                     *
 * l_v               --> length of the interaction list                      *
 * event[end+1:l_v]  --> interaction list starting from absorption           *
 *****************************************************************************/

#include "stdio.h"
#include "binio.h"
#include "event.h"

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

static int i,j,k;
static int v_l;
static int p_l;

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//picks an event from the data with index ndx
int pick_event(double event[], double data[], int n_proc[], 
	       int n_ev, int ndx, int verb){
  v_l = 7; p_l = 13;
  int f_ndx=0,l_ndx=0;
  int list[15];
  double time[15];
  int added;
  
  if(ndx<n_ev){
    if(verb) print_one(data,n_proc,ndx);
    for(i=0;i<ndx;i++) f_ndx += v_l + n_proc[i]*p_l;
    event[0] = (double) n_proc[ndx];
    for(i=0;i<event[0];i++){
      event[1+i*4] = data[f_ndx+v_l+i*p_l+12];
      for(j=0;j<3;j++){
	event[2+i*4+j] = data[f_ndx+v_l+i*p_l+3+j];
      }
      if(i==0){
	time[1]=data[f_ndx+v_l+i*p_l+2];
	list[1]=0;
      }else{
	time[0]=data[f_ndx+v_l+i*p_l+2];
	list[0]=i;
	added=0;
	for(j=1;j<i+1;j++){
	  if(time[0]>time[j]){
	    added=1;
	    for(k=i;k>=j;k--){
	      time[k+1]=time[k];
	      list[k+1]=list[k];
	    }
	    time[j]=time[0];
	    list[j]=list[0];
	    break;
	  }
	}
	if(!added){
	  list[i+1] = i;
	  time[i+1] = data[f_ndx+v_l+i*p_l+2];
	}
      }
    }
    for(j=0;j<event[0];j++) event[1+(int)event[0]*4+j] = (double)list[j+1];
    if(verb){
      printf("Interaction order :\n");
      printf("* ");
      for(j=0;j<event[0];j++){
	printf(" %d -->",(int)event[1+(int)event[0]*4+j]);
      }
      printf(" Sky\n");
    }   
    return 0;
  } else {
    printf("\n\nERROR!!\n-->Choose a valid index...\n");
    return 1;
  }

}

int event_order(int list[], double event[], double data[]){

  return 0;
}

/*****************************************************************************/
