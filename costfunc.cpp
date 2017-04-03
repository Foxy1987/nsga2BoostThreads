#include <cstdlib>
#include <math.h>
#include <iostream>
#include "globals.h"


void costfunc(double *x, double *obj, double *constr) {
  obj[0] = pow(((x[0] - Z0)/Z0), 2);
  obj[1] = pow(((x[1] - FMAX)/FMAX), 2);
  obj[2] = pow(((x[2] - ZMAX)/ZMAX), 2);
  obj[3] = pow(((x[3] - Z10)/Z10), 2);
  obj[4] = pow(((x[4] - FWIDTH1)/FWIDTH1), 2);
  obj[5] = pow(((x[5] - FWIDTH2)/FWIDTH2), 2);
  obj[6] = pow(((x[6] - PHIZERO)/PHIZERO), 2);
  obj[7] = pow(((x[7] - PHI0PT1)/PHI0PT1), 2);
  obj[8] = pow(((x[8] - PHI2)/PHI2), 2);
  obj[9] = pow(((x[9] - PHI4)/PHI4), 2);
  /*
    //  Z0 
    double a, b;
    //a = 8.6035-x[0]; b = x[0]-7.7841;
    a=8.6-x[0]; b = x[0]-7.8;
    //a=9-x[0]; b = x[0]-8;

    if (a < 0) {
        constr[0] = a;
    } else if (b < 0) {
        constr[0] = b;
    } else {
        constr[0] = 0;
    }

    // FMAX 
    a = 1.1-x[1]; b = x[1]-1;
    //a = 1.1-x[1]; b = x[1]-0.9;
    if (a < 0) {
        constr[1] = a;
    } else if (b < 0) {
        constr[1] = b;
    } else {
        constr[1] = 0;
    }

    // ZMAX 
    //a = 5.7924-x[2]; b = x[2] - 5.2407;
    a = 13.8-x[2]; b = x[2]-13.5; 
    if (a < 0) {
        constr[2] = a;
    }else if (b < 0) {
        constr[2] = b;
    } else {
        constr[2] = 0;
    }

    // Z10
    a = 10-x[3]; b=x[3]-9;
    if (a < 0) {
        constr[3] = a;
    }else if (b < 0) {
        constr[3] = b;
    } else {
        constr[3] = 0;
    }
    

    // Fwidth1 
    a = 0.42-x[4]; b=x[4]-0.38;
    if (a < 0) {
        constr[4] = a;
    }else if (b < 0) {
        constr[4] = b;
    } else {
        constr[4] = 0;
    }

    // Fwidth2
    a = 2.6-x[5]; b=x[5]-2.3;
    if (a < 0) {
        constr[5] = a;
    }else if (b < 0) {
        constr[5] = b;
    } else {
        constr[5] = 0;
    }

    //  phasonance
    a = 1.1-x[6]; b =x[6]-1;
    if (a < 0) {
	constr[6] = a;
    }
    else if (b < 0) {
	constr[6] = b;
    }
    else {
	constr[6] = 0;
    }

    //a = 0.15-x[7]; b=x[7]-0.1367;
    a = 0.2-x[7]; b=x[7]-0;
    if (a < 0) {
       constr[7] = a;
    } 
    else if (b < 0) {
	constr[7] = b;
    }
    else {
    	constr[7] = 0;
    }

    a = -0.1-x[8]; b=x[8]+0.4;
    if (a < 0) {
	constr[8]=a;
    }
    else if (b < 0) {
        constr[8] = b;
    }
    else {
        constr[8] = 0;
    }

    a = -0.3-x[9]; b=x[9]+0.45;
    if (a < 0) {
        constr[9]=a;
    }
    else if (b < 0) {
	constr[9] = b;
    }
    else {
        constr[9] = 0;
    }
*/

  
    //  Z0 
    double a, b;
    //a = 8.6035-x[0]; b = x[0]-7.7841;
    a=8.6-x[0]; b = x[0]-7.8;
    //a=9-x[0]; b = x[0]-8;

    if (a < 0) {
        constr[0] = a;
    } else if (b < 0) {
        constr[0] = b;
    } else {
        constr[0] = 0;
    }

    // FMAX 
    a = 1.1-x[1]; b = x[1]-1;
    //a = 1.1-x[1]; b = x[1]-0.9;
    if (a < 0) {
        constr[1] = a;
    } else if (b < 0) {
        constr[1] = b;
    } else {
        constr[1] = 0;
    }

    // ZMAX 
    //a = 5.7924-x[2]; b = x[2] - 5.2407;
    a = 13.6-x[2]; b = x[2]-13.1; 
    if (a < 0) {
        constr[2] = a;
    }else if (b < 0) {
        constr[2] = b;
    } else {
        constr[2] = 0;
    }

    // Z10
    a = 10.5-x[3]; b=x[3]-9.5;
    if (a < 0) {
        constr[3] = a;
    }else if (b < 0) {
        constr[3] = b;
    } else {
        constr[3] = 0;
    }
    

    // Fwidth1 
    a = 0.5-x[4]; b=x[4]-0.38;
    if (a < 0) {
        constr[4] = a;
    }else if (b < 0) {
        constr[4] = b;
    } else {
        constr[4] = 0;
    }

    // Fwidth2
    a = 2.6-x[5]; b=x[5]-2.4;
    if (a < 0) {
        constr[5] = a;
    }else if (b < 0) {
        constr[5] = b;
    } else {
        constr[5] = 0;
    }

    //  phasonance
    a = 1.1-x[6]; b =x[6]-1;
    if (a < 0) {
	constr[6] = a;
    }
    else if (b < 0) {
	constr[6] = b;
    }
    else {
	constr[6] = 0;
    }

    //a = 0.15-x[7]; b=x[7]-0.1367;
    a = 0.2-x[7]; b=x[7]-0;
    if (a < 0) {
       constr[7] = a;
    } 
    else if (b < 0) {
	constr[7] = b;
    }
    else {
    	constr[7] = 0;
    }

    a = -0.25-x[8]; b=x[8]+0.4;
    if (a < 0) {
	constr[8]=a;
    }
    else if (b < 0) {
        constr[8] = b;
    }
    else {
        constr[8] = 0;
    }

    a = -0.3-x[9]; b=x[9]+0.45;
    if (a < 0) {
        constr[9]=a;
    }
    else if (b < 0) {
	constr[9] = b;
    }
    else {
        constr[9] = 0;
    }

  }

