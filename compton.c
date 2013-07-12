/*****************************************************************************\
 * Yigit Dallilar 10.07.2013                                                 *
 * DTU-Space                                                                 *
 * Compton telescope looks for interaction order                             *
\*****************************************************************************/

#include "math.h"
#include "compton.h"

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

float rest_e = 511;

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

double compt_angle(double post, double pre){
  
  return acos(rest_e*(pre-post)/(pre*post)-1);
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int compt_match(double match[], double path[], double event[], int n_ev,
		int verb){ 


}
/*****************************************************************************/
