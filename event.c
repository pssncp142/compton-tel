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
 * event[end+l_v+1]  --> if the event is complete or not                     *
 *****************************************************************************/

#include "stdio.h"
#include "binio.h"
#include "event.h"
#include "math.h"

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

static int i,j,k;
static int v_l=7;
static int p_l=13;

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//picks an event from the data with index ndx
int pick_event(double event[], double data[], int n_proc[], 
	       int n_ev, int ndx, int verb){
  //v_l = 7; p_l = 13;
  int f_ndx=0,l_ndx=0;
  double tot_en=0;
  int n_p; 
  
  if(ndx<n_ev){
    if(verb) print_one(data,n_proc,ndx);
    for(i=0;i<ndx;i++) f_ndx += v_l + n_proc[i]*p_l;
    n_p = n_proc[ndx];
    event[0] = (double)n_p;
    for(i=0;i<n_p;i++){
      tot_en += data[f_ndx+v_l+i*p_l+12];
      event[1+i*4] = data[f_ndx+v_l+i*p_l+12];
      for(j=0;j<3;j++){
	event[2+i*4+j] = data[f_ndx+v_l+i*p_l+3+j];
      }
    }
    if(tot_en>data[f_ndx+6]-0.001){ 
      event[1+n_p*5]=1;
    } else {
      event[1+n_p*5]=0;
    }

    event_order(event,data,n_proc,ndx,verb);
    return 0;
  } else {
    printf("\n\nERROR!!\n-->Choose a valid index...\n");
    return 1;
  }

}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/
//adds event order to event array just after the events are finished.

int event_order(double event[], double data[], int n_proc[],
		int ndx, int verb){

  double time[15];
  int list[15];
  int added;
  int f_ndx=0;
  int n_p = (int)floor(event[0]); 

  for(i=0;i<ndx;i++) f_ndx += v_l + n_proc[i]*p_l;

  for(i=0;i<n_p;i++){
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

  for(j=0;j<n_p;j++) event[1+n_p*4+j] = (double)list[j+1];
  if(verb){
    printf("Interaction order :\n");
    printf("* ");
    if(event[1+n_p*5]==0) printf(" Out -->");
    for(j=0;j<n_p;j++){
      printf(" %d -->",(int)event[1+n_p*4+j]);
    }
    printf(" Sky\n");
  }   

  return 0;
}

/*****************************************************************************/
