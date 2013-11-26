#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>

/*----------------------------------------------------------*/
/* time in sec                                              */
/*----------------------------------------------------------*/
double util_gettime(){
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)tp.tv_sec+(0.000001)*tp.tv_usec);
}
