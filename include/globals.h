#ifndef __GLOBALS_H
#define __GLOBALS_H

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>

/*popsize must be a multiple of 4 */
//#define POP_SIZE 100
#define POP_SIZE 200
#define NUM_PARS 8 
#define NUM_OBJS 10 
#define NUM_GEN 2000
#define NUM_CON 10
//#define NUM_CON 9 


/*#define Z0 8.1938 
#define FMAX 1.0341
#define ZMAX 13.71
#define Z10 9.5628
#define FWIDTH1 0.4083
#define FWIDTH2 2.4984
#define PHIZERO 1.05
#define PHI2 -0.41
#define PHI0PT1 0.1
#define PHI4 -0.38
*/
#define Z0 8.15
#define FMAX 1.03
#define ZMAX 13.25
#define Z10 10
#define FWIDTH1 0.43
#define FWIDTH2 2.57
#define PHIZERO 1.06
#define PHI2 -0.37
#define PHI0PT1 0.11
#define PHI4 -0.37

#define INF 1.0e14                              
#define EPS 1.0e-14
/* crossover distribution index between 5 and 20 */
// #define ETA_C 20
/* mutation distribution index between 5 and 50 */
// #define ETA_M 20 
// #define P_CROSS 0.7 
// #define P_MUT 0.7 

#define ETA_C 5 
#define ETA_C 5
// mutation distribution index between 5 and 50 
#define ETA_M 20
#define P_CROSS 0.8
#define P_MUT 0.7

#define NUM_THREADS 10 
#define NUM_SIMS 20 
#define PI 3.141592653589793

void costfunc(double *, double *, double *);

#endif

