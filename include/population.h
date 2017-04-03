#ifndef __POPULATION_H
#define __POPULATION_H

#include <vector>
#include <iterator>
#include <algorithm>
#include <math.h>
#include <map>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include "individual.h"
#include "worker.h"

using namespace std;
using namespace boost;

class Population {
  public:
    static double *min_var;
    static double *max_var;
    Individual **ind;
    int _popsize;
    int nreal;    /* number of real paramaters */
    int nobj; /* number of objective functions */
    int ncon; /* number of constraints */
    Random *rand; /* random number generator */
    int numFronts;
    
    boost::mutex lock;
   
    static int ncross;
    vector<vector<Individual *> > fronts;
   
    void simulate();
        /* inner class whose instances are intended to be executed by a thread
         * Each thread will run a number of simulations
         */
  
    Population(int popsize, int nreal, int nobj, int ncon);
    ~Population();
    void initialize();
    void print_pop(ostream& os);
    void printConstraintVals(ostream& os);
    void readLimits(ifstream &ifs);
    void evaluate();
    void nondominated_sort();
    void  assignCrowdingDistance(vector<Individual *> &, int frontsize);
    void generateNewParent(Population *);
    void quicksortOnObjective(vector<Individual *> &front, int objCount, int left, int right);
    int randPartition(vector<Individual *> &, int objCount, int left, int right);
    void swap(vector<Individual *> &, int a, int b );
    void quicksortOnCrowdingDistance(vector<Individual *> &, int left, int right);
    int randPartitionDist(vector<Individual *>&, int left, int right);
    void select (Population *new_pop);
    Individual* tournament(Individual *ind1, Individual *ind2);
    void sbxCrossover (Individual *parent1, Individual *parent2, Individual *child1, Individual *child2);
    void mutate_pop();
    void merge(Population *, Population *);
    void printObjectives(ostream& os);
    void printAttributes(ostream& os);
    void qsortOnRank(int left, int right);
    int randPartitionRank(int left, int right);
    int getPopSize();
    void printGoodSolutions(ostream& os, int genNum);
};
#endif
