/*****************************************************************************\
 * Yigit Dallilar 12.07.2013                                                 *
 * DTU-Space                                                                 *
 * Compton telescope library for events                                      *
\*****************************************************************************/

#include "stdio.h"
#include "binio.h"
#include "event.h"

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int i,j,k;
int v_l;
int p_l;

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int pick_event(double event[], double data[], int n_proc[], 
	       int n_ev, int ndx, int verb){
  v_l = 7; p_l = 13;
  int f_ndx=0,l_ndx=0;
  
  if(ndx<n_ev){
    if(verb) print_one(data,n_proc,ndx);
    for(i=0;i<ndx;i++) f_ndx += v_l + n_proc[i]*p_l;
    event[0] = (double) n_proc[ndx];
    for(i=0;i<event[0];i++){
      event[1+i*4] = data[f_ndx+v_l+i*p_l+12];
      for(j=0;j<3;j++){
	event[2+i*4+j] = data[f_ndx+v_l+i*p_l+3+j];
      }
    }
    return 0;
  } else {
    printf("\n\nERROR!!\n-->Choose a valid index...\n");
    return 1;
  }

}

/*****************************************************************************/
