/*****************************************************************************\
 * Yigit Dallilar 10.07.2013                                                 *
 * DTU-Space                                                                 *
 * Compton telescope data in out library                                     *
\*****************************************************************************/

#include "stdio.h"
#include "binio.h"

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int read_bin(char fname[], double data[], int n_proc[]){

  int i,j,ndx,tmp,n_ev;
  char a;
  FILE* f;
  f = fopen(fname,"r");

  ndx = 0;
  n_ev = 0;
  while(getc(f) != EOF){
    fread(&tmp,sizeof(int),1,f);
    n_proc[n_ev] = tmp;
    for(i=0;i<7;i++) fread(&data[ndx++],sizeof(double),1,f);
    for(i=0;i<n_proc[n_ev];i++) {
      for(j=0;j<2;j++){
	fread(&tmp,sizeof(int),1,f);
	data[ndx++] = (double) tmp;
      }
      for(j=0;j<11;j++) fread(&data[ndx++],sizeof(double),1,f);
    }
    n_ev++;
  }
  fclose(f);
  
  return n_ev;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int print_one(double data[], int n_proc[], int ndx){

  int i,j;
  int st_ndx=0;

  for(i=0;i<ndx;i++) st_ndx += 7 + n_proc[i]*13;

  printf("\n\nEvent Index    : %2d\n",ndx);
  if(n_proc[ndx]==0){
    printf("Process Count  : %2d\n",0);
    printf("***No process to display.***\n");
  } else {
    printf("Process Count  : %2d\n",n_proc[ndx]);
    printf("Vertex Info    :\n");
    printf("  -Position  : %7.3f , %7.3f , %7.3f\n",data[st_ndx],data[st_ndx+1],data[st_ndx+2]);
    printf("  -Momentum  : %7.3f , %7.3f , %7.3f\n",data[st_ndx+3],data[st_ndx+4],data[st_ndx+5]);
    printf("  -Kinetic   : %7.3f\n\n",data[st_ndx+6]);
    st_ndx += 7;
    printf("Process Info   :\n");
    printf(" No Proc Det Time(s)         Position(mm)            PreMomentumDir           PostMomentumDir        EnergyDeposit(keV)\n");
    for(i=0;i<n_proc[ndx];i++){
      printf("*%2d",i);
      if(data[st_ndx]==0){
	printf(" %s","cmpt");
      } else if(data[st_ndx]==1){
	printf(" %s","phot");
      } else if(data[st_ndx]==2){
	printf(" %s","rayl");
      }
      if(data[st_ndx+1]==0){
	printf(" %s","ins");
      }else if(data[st_ndx+1]==1) {
	printf(" %s","out");
      }
      printf(" %6.3f  (%7.3f,%7.3f,%7.3f) (%7.3f,%7.3f,%7.3f) (%7.3f,%7.3f,%7.3f) %15.3f \n",
	     data[st_ndx+2],data[st_ndx+3],data[st_ndx+4],data[st_ndx+5],data[st_ndx+6],data[st_ndx+7],data[st_ndx+8],
	     data[st_ndx+9],data[st_ndx+10],data[st_ndx+11],data[st_ndx+12]);
      st_ndx += 13;
    }
  }

  return 0;
}

/*ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo*/

int print_all(double data[], int n_proc[], int n_ev){
  
  int i;
  for(i=0;i<n_ev;i++){
    print_one(data,n_proc,i);
  }

  return 0;
}

/*****************************************************************************/
