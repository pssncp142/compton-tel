/*****************************************************************************\
 * Yigit Dallilar 10.07.2013                                                 *
 * DTU-Space                                                                 *
 * Compton telescope - returns statistics about the data                     *
\*****************************************************************************/

#include "stdio.h"
#include "stat.h"

static int f_ndx,l_ndx;
static int i,j,k;
static int fnd;
static int v_l;
static int p_l;

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int stat_complete(double stat[], double data[], int n_proc[], int n_ev, int verb){

  v_l =  7;
  p_l = 13;
  f_ndx = 0;
  
  double ener;

  stat[0] = 0;
  stat[1] = 0;
  stat[2] = 0;

  for(i=0;i<n_ev;i++){
    fnd = 0;
    ener = 0;
    l_ndx = f_ndx + v_l + n_proc[i]*p_l;
    for(j=0;j<n_proc[i];j++){
      ener += data[f_ndx+v_l+p_l*j+12];
    }
    if(data[f_ndx+6] != ener) fnd = 1;
    if(ener == 0) fnd = 2;
    stat[fnd]++;
    f_ndx = l_ndx;
  }

  stat[0] /= n_ev;
  stat[1] /= n_ev;
  stat[2] /= n_ev;

  printf("From %3d of events...\n",n_ev);
  printf("All energy detected events : %4.2f %%\n",stat[0]*100);
  printf("Partial detected events    : %4.2f %%\n",stat[1]*100);
  printf("No interaction             : %4.2f %%\n",stat[2]*100);
  
  
  return 0;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int stat_nofproc(double stat[], double data[], int n_proc[], int n_ev, int verb){

  v_l =  7;
  p_l = 13;

  return 0;
}

/*****************************************************************************/
