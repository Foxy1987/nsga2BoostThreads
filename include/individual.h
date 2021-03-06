#ifndef __INDIVIDUAL_H
#define __INDIVIDUAL_H

#include <iostream>
//#include <cstdlib>
#include "random.h"

#include <stdio.h>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>    // std::copy
#include <vector>
#include "globals.h"

#include <boost/scoped_array.hpp>
#include <boost/lexical_cast.hpp>
using namespace std;
using namespace boost;

class Individual {
  public:
    double *objs; /* fitness functions */
    double *vars; // stores params
    double *xreal; //stores zfstats
    double *constr;
    double crowding_distance;
    int rank;
    int numvars;
    int nobj;
    int ncon;
    int ndom;    // domination count
    vector<Individual *> ss; // set of solutions dominated by this solution
    
    double constr_violation;
    int index;
    
    Individual(int nreal, int nobj, int ncon);
    ~Individual();
    void initialize(vector<string> &v);
    void initialize(Random *rand, double *, double *);
    void evaluate();
    int checkDominance(Individual *);
    void mutate(Random *rand, double *, double *);
    void simulate(int id);
    void copy(Individual *);
};
#endif



