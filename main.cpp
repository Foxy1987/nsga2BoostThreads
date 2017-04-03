#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cstring>
#include "population.h"
#include <iomanip>

double* Population::min_var = new double[NUM_PARS];
double* Population::max_var = new double[NUM_PARS];

int main(int argc, char *argv[]) {
    
    Population *parent_pop = new Population(POP_SIZE, NUM_PARS, NUM_OBJS, NUM_CON);
    Population *child_pop = new Population(POP_SIZE, NUM_PARS, NUM_OBJS, NUM_CON);
    Population *mixed_pop = new Population(2*POP_SIZE, NUM_PARS, NUM_OBJS, NUM_CON);
    ifstream limitsfile("limits2.dat");
  
     parent_pop->readLimits(limitsfile);
     
     limitsfile.close();
        /* initialize parent population with randomly generated
         * numbers using boost libraries */
    parent_pop->initialize();
        /* send initial population to a file */
        /* initial parallel simulations */
    parent_pop->simulate();
    cout << "Finished Simulating\n";   
        // evaluate entire population
    parent_pop->evaluate();
    
    int flag = 0;
    /* perform non-dominated sorting and assign crowding distance */
    parent_pop->nondominated_sort(); 
   
    ofstream os1, os2, os3;
    setprecision(5);

    os1.open("final_pop.dat", ios::app);
    parent_pop->print_pop(os1);

    //os2.open("/afs/cad/u/d/m/dmf6/HD/gorgonscratch/BestSolutions.dat", ios::app);
    os2.open("BestSolutions.dat", ios::app);
       /*  Start evolution */
     for (int i = 1; i < NUM_GEN; i++) {
       cout << "GEN NUM " << i << endl;
       cout << "Starting generation " << i << "\n";

       /* select parents and produce offspring by SBX crossover */
       parent_pop->select(child_pop);
       /*mutate children */
       child_pop->mutate_pop();
       /* this is done in parallel */
       child_pop->simulate();
       
       /* evaluate how good of a match to target zf profile */
       child_pop->evaluate();
             
       /* populate mixed_pop with parent and child populations */
       mixed_pop->merge(parent_pop, child_pop);       
       
       mixed_pop->nondominated_sort();  //non-dominated sort on all fronts
  
       mixed_pop->generateNewParent(parent_pop);
       //cout << "FINISHED GENERATING NEW PARENT\n";
       
       /* print the good solutions that were found in this generation */
       parent_pop->printGoodSolutions(os2, i);

       /* print parameters */
       cout << "Generations finished, now sending solutions to final_pop.dat" << endl;
       parent_pop->print_pop(os1); 
     }

     os1.close();
     os2.close();
     return 0;
}
