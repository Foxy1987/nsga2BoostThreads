#include "population.h"
#include <algorithm>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

int Population::ncross = 0;

 bool compareDistance(const Individual *a, const Individual *b) {
    return (a->crowding_distance > b->crowding_distance);
}

/* inline constructor */
Population::Population(int popsize, int nreal, int nobj, int ncon) : _popsize(popsize) {
        /* need to seed from the random number generator for the entire population */
    rand = new Random();
    rand->setSeed(static_cast<unsigned int>(time(0)));
    this->nreal = nreal;
    this->nobj = nobj;
    ind = new Individual*[_popsize];
  
    for (int i = 0; i < _popsize; i++) {
        ind[i] = new Individual(nreal, nobj, ncon);
        ind[i]->index = i;
    }
    ncross = 0; 
    this->ncon = ncon;
    numFronts = 0;
}

Population::~Population() {
    for (int j = 0; j < _popsize; j++) {
        delete ind[j];
    }
    delete [] ind;
    delete rand;
}

/* initialize all of the individuals already created in constructor
 * with random parameters */ 
void Population::initialize() {
    vector<string> tokenList;
    string line; 
    int i = 0;
    char_separator<char> sep("\t");
    int flag = 0; 
    if (flag) {
        ifstream init("init.dat");

        if (init.is_open()) {     
            while (i < _popsize) {
	            getline(init, line);
	            cout << line << "\n\n\n\n\n";
	            tokenizer<char_separator<char> > tokens(line, sep);
	            BOOST_FOREACH(string t, tokens) {
	                tokenList.push_back(t);
//	                cout << t << "." << "\n\n\n\n\n\n";
	            }
                ind[i]->initialize(tokenList);
                tokenList.clear();
//	            //split(tokenList, line, is_any_of("\t"), token_compress_on);
	            i++;
            }
        }
        init.close();  
    }
    else {
        for (int i = 0; i < _popsize; i++) {
            //cout << "IND i\n";
            ind[i]->initialize(rand, min_var, max_var);     
            cout  <<  i  <<  "\n";
        }
    }
    cout  <<  "Finished Initializing\n";
}

void Population::print_pop(ostream& os) {
        //os  << "Population size = " << _popsize << "\t" << " Number of Parameters being optimized = " << nreal << endl;
 
    for (int i = 0; i < _popsize; i++) {
        for (int j = 0; j < nreal; j++) {
            os << " " << ind[i]->vars[j] << " ";
        }
        for(int j = 0; j < nobj; j++) {
            os << " " << ind[i]->xreal[j] << " ";
        }
        os << "\n";   
    }
    os << endl;
    
}

void Population::printObjectives(ostream& os) {
    for (int i = 0; i < _popsize; i++) {
        for(int j = 0; j < nobj; j++) {
            os << " " << ind[i]->objs[j] << " ";
        }
        os << "\n";   
    }
    os << endl;
}

void Population::printAttributes(ostream& os) {
    for (int i = 0; i < _popsize; i++) {
        for(int j = 0; j < nobj; j++) {
            os << " " << ind[i]->xreal[j] << " ";
        }
        os << "\n";   
    }
    os << endl;
}



void Population::printGoodSolutions(ostream& os, int genNum) {

	int flag = 1; int count = 0;
	for (int i = 0; i < _popsize; i++) {
		if (ind[i]->constr_violation < 0) {
			flag = 0;
		}
		if (flag) {
			for (int k = 0; k < nreal; k++) {
				os << ind[i]->vars[k] << "\t";
			}
			for (int l = 0; l < nobj; l++) {
				os << ind[i]->xreal[l] << "\t";
			}
			os << genNum << endl;
			count++;
		}
	}
	if (count == 0) {
		os << "No Solutions were found at generation number: " << genNum << endl;
	}
}

void Population::evaluate() {
    for (int j = 0; j < _popsize; j++) {
            // each individual in the popultion will be evaluated against each objective
        Individual *current_ind = ind[j];
        current_ind->evaluate();
    }
}




void Population::nondominated_sort() {
        //clear fronts vector for the current generation
        cout << "SORTING!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    fronts.clear();
    vector< Individual *> front;    /* first front */
  
    int flag = 0;
    
    for (int p = 0; p < _popsize; p++) {
        Individual *ind_p = ind[p];
        ind_p->ss.clear();
        ind_p->ndom=0;
        for (int q = 0; q < _popsize; q++) {
            if (p != q) {
                Individual *ind_q = ind[q];
                // check to see if p dominates q
                flag = ind_p->checkDominance(ind_q);      
//                cout  <<  "FLAG IS "  <<  flag  <<  "\n";
                if (flag==1) {
                    ind_p->ss.push_back(ind_q);
                } else if (flag == -1){
                //update domination count for individual 
                    ind_p->ndom++;
                }
            }
        }
        if (ind_p->ndom == 0) {
                /* If there are no individuals that dominate ind_p then it belongs to the
                 *  first front. Each individual has a member variable called rank
                 */
                //cout << "No individuals dominate individual " << ind_p->index << endl;
            
            ind_p->rank = 1;
                //store individual in first front because it is 'non-dominated' by all solutions
            front.push_back(ind_p);
                /* store index of element too, first front element indices come first obviously */
            cout  <<  front.size()  <<  "\n";
        }
    }
    /* store first front */
    fronts.push_back(front);
    assignCrowdingDistance(front, front.size());

//    int start = front.size();
    
    numFronts = 2; //representing the first front
    vector<Individual *> q; /* used to store members of the next front */
//    cout << "HELLO!\n";
    while(!front.empty()) {
            //cout << "Creating a new front..." << endl;
        q.clear();
            /* iterate through all members of the first front */
        vector<Individual *>::iterator it = front.begin();
        //vector.end returns an iterator to the element past the end of the sequence
        while (it !=front.end()) {
            Individual *indi = *it;
            vector<Individual *> s = indi->ss;
             /* iteratre through the set of solutions dominated by it 
    		 s is the set of individuals dominated by ind*/
            foreach(Individual *qi, s) {
                qi->ndom -=1;
                // if ndom = 0, then none of the individuals in the subsequent fronts would dominate q 
                if(qi->ndom == 0) {
                        //cout << "domination count for individual with idx " << qi->index  << " is " << qi->ndom << endl;
                    qi->rank = numFronts;
                    q.push_back(qi);
                        //cout << "The number of individuals in q is " << q.size() << "\n";
                        /* store index of individual whose domination
                           count has reached 0. This individual will
                           not be visited again */
                }
            }
	    it = front.erase(it);
        }
        
        if (!q.empty()) {
	    fronts.push_back(q);
            assignCrowdingDistance(q, q.size());
            for (int i =0; i < q.size(); i++) {
		cout << "crowding distance in q[i] is " << q[i]->crowding_distance << "\n";
	    }
	    numFronts++; //increment i to identify the next front
        }
        front.swap(q);
    }
    qsortOnRank(0, _popsize-1);
    for (int i = 0; i < _popsize; i++) {
         cout  <<  ind[i]->rank  << "\t" << ind[i]->crowding_distance << "\n";
    }
    cout  <<  "\n";
    cout << "FINISHED NON-DOMINATED SORTING WITH " << fronts.size() << "sets\n";
   
}

/*randomized quicksort to sort front based on objection objcount */
void Population::qsortOnRank(int left, int right) {
        //initialize pivot index
    int pivot;
    if (left < right) {
        pivot = randPartitionRank(left, right);
        qsortOnRank(left, pivot-1);
        qsortOnRank(pivot+1, right);
    }
}

/* we don't pass the population ind because it is a member of this class */
int Population::randPartitionRank(int left, int right) {
    int pivot;
    int randIdx = rand->nextInt(left, right);

    std::swap(ind[right], ind[randIdx]);

    pivot = ind[right]->rank;
    
     int k = left - 1;
     for (int j = left; j < right; j++) {
         if(ind[j]->rank <=pivot) {
              k++;
              std::swap(ind[k], ind[j]);
         }
     }
     std::swap(ind[k+1], ind[right]);
    return k+1;
}


/* notice that the population members aren't actually being switched
 * just the indices array and this is used to index the population
 * when calculating the crowding distance */
void Population::assignCrowdingDistance(vector<Individual *> &front, int frontsize){
  /* reset crowding distance for all solutions */
  std::vector<Individual *>::const_iterator it;
  Individual *i;
  for (it = front.begin(); it != front.end(); ++it) {
      i = *it;
      i->crowding_distance = 0;
  }
  if (front.size() == 1){
    cout << front.front()->crowding_distance << "\n";
    cout << "The front size is 1\n";
    front.front()->crowding_distance = INF; 
    return;
  }
  if (front.size() == 2) {
    cout << "THE FRONT SIZE IS 2\n";  
    front.front()->crowding_distance = INF; 
    front.back()->crowding_distance = INF;
    return;
  }

        /* for each objective call randqsort to sort the set in worse
         * order of current obejctive (return sorted indices array) */
    // now iterate through vector of Individual *    

  cout << "THE FRONT SIZE IS " << frontsize << "\n";
  for (int i = 0; i < nobj; i++) {
      
//            //using randomized quicksort, sort the front indices array
//            // in order of objective i use the sorted indices array to
//            // index the population when calculating crowding distance
//            // for all members in front
//  for each objective, the individuals are sorted in ascending order based on their value for this objective
    
    quicksortOnObjective(front, i, 0, frontsize-1);
    
//
//    /* assign infinite distance to boundary solutions i.e., front and back */
    front.front()->crowding_distance = INF;
    front.back()->crowding_distance = INF;
//    /* maximum objective value */
    double fmax = front[frontsize-1]->objs[i];
    
   
//    /* minimum objective value */
    double fmin = front[0]->objs[i];
         
//    /* all other solutions 2-(n-1) */
    double nextObj, previousObj;
    for (int j = 1; j < frontsize-1; j++) {
      nextObj = front[j+1]->objs[i];
      previousObj =  front[j-1]->objs[i];
      
      if ((fmax - fmin) == 0) {
          front[j]->crowding_distance = INF;
      }
      else if ((nextObj - previousObj) == 0) {
          front[j]->crowding_distance = INF;
      }
      else {
            // what if nextObj and previousObj are the same?
	    cout << front[j]->xreal[2] << "\n";
	    cout << previousObj << " " << nextObj << " " << fmin << " " << fmax << " " << (nextObj - previousObj)/( fmax - fmin) << "\n";
	    front[j]->crowding_distance += (nextObj - previousObj)/( fmax - fmin);
	    cout << "CROWDING DISTANCE 1: " << front[j]->crowding_distance << "\n";
      }    
    }
  }
  
  for (int k = 0; k < frontsize; k++) {
    if (front[k]->crowding_distance != INF) {
      front[k]->crowding_distance = (front[k]->crowding_distance)/nobj;
    }
    cout << "CROWDING DISTANCE 2 " << front[k]->crowding_distance << "\n";
  }
  cout << endl;
}

/*randomized quicksort to sort front based on objective objcount */
void Population::quicksortOnObjective(vector<Individual *> &front, int objCount, int left, int right) {
        //initialize pivot index
    int pivot;
    if (left < right) {
            //pass objective number
        pivot = randPartition(front, objCount, left, right);
        quicksortOnObjective(front, objCount, left, pivot-1);
        quicksortOnObjective(front, objCount, pivot+1, right);
    }
}

/* we don't pass the population ind because it is a member of this class */
int Population::randPartition(vector<Individual *> &front, int objCount, int left, int right) {
    double pivot;
    int randIdx = rand->nextInt(left, right);
    swap(front, right, randIdx);
        /* use random index as pivot and use indices to index the
         * population and objCount to index objective array*/
    pivot = front[right]->objs[objCount];
    int k = left - 1;
    for (int j = left; j < right; j++) {
        if(front[j]->objs[objCount] <=pivot) {
            k++;
            swap(front, k, j);
        }
    }
    swap(front, k+1, right);
    return k+1;
}

void Population::swap(vector<Individual *> &array, int a, int b ) {
    Individual *temp = array[ a ];
    array[ a ] = array[ b ];
    array[ b ] = temp;
}

bool equals (Individual *i, Individual *j) {
  return (i->index == j->index);
}

void Population::select (Population *new_pop) {
 //cout << "STARTING SELECTION\n";  
  int c1, c2; 
  vector<Individual *> matingPool(_popsize);
  Individual *parent;
  for (int i = 0; i < _popsize; i++) {
  	c1=rand->nextInt(0, _popsize-1);
  	c2=rand->nextInt(0, _popsize-1);
  	while (c1 == c2) {
		c2 = rand->nextInt(0, _popsize-1);
  	}
//	cout << "C1 = " << c1 << " C2 = " << c2 << "\n";
	parent = tournament(ind[c1], ind[c2]);
	matingPool[i] = parent;
  }
// cout << " FINISHED BUILDING MATING POOL!!!!!!!!!!!!!!\n";
  	// now loop through mating pool to generate child population
  int k = 0;
  while (k < _popsize-1) {
    sbxCrossover (matingPool[k], matingPool[k+1], new_pop->ind[k], new_pop->ind[k+1]);
   k = k + 2;
  }


}

Individual* Population::tournament (Individual *ind1, Individual *ind2) {
    double crowd1 = ind1->crowding_distance;
    double crowd2 = ind2->crowding_distance;
    
    int rank1 = ind1->rank;
    int rank2 = ind2->rank;

    if (rank1 < rank2){
        return (ind1);
    }
    else if (rank2 < rank1) {
        return (ind2);
    }
    else if (crowd1 > crowd2) {
        return(ind1);
    }
    else if (crowd2 > crowd1) {
        return(ind2);
    }
    else if (rand->nextDouble(0, 1) <= 0.5) {
        return(ind1);
    }
    else {
        return(ind2);
    }
}

/* Routine for real variable SBX crossover */
void Population::sbxCrossover (Individual *parent1, Individual *parent2, Individual *child1, Individual *child2) {
    double random;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;
    
    if (rand->nextDouble(0, 1)<=P_CROSS) {
        for (int i=0; i<NUM_PARS; i++) {
            if (rand->nextDouble(0, 1) <= 0.5) {
                 if (fabs(parent1->vars[i]-parent2->vars[i]) > EPS) {
                     if (parent1->vars[i] < parent2->vars[i]) {
                         y1 = parent1->vars[i];
                         y2 = parent2->vars[i];
                     }
                     else {
                         y1 = parent2->vars[i];
                         y2 = parent1->vars[i];
                     }
                     yl = min_var[i];
                     yu = max_var[i];

                     random = rand->nextDouble(0, 1);
                     beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                     alpha = 2.0 - pow(beta,-(ETA_C+1.0));
                     if (random<= (1.0/alpha)) {
                         betaq = pow ((random*alpha),(1.0/(ETA_C+1.0)));
                     }
                     else {
                         betaq = pow ((1.0/(2.0 - random*alpha)),(1.0/(ETA_C+1.0))); 
                     }
                     c1 = 0.5*((y1+y2)-betaq*(y2-y1));                
                     beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                     alpha = 2.0 - pow(beta,-(ETA_C+1.0));
                     if (random <= (1.0/alpha)) {
                         betaq = pow ((random*alpha),(1.0/(ETA_C+1.0)));
                     }
                     else {
                         betaq = pow ((1.0/(2.0 - random*alpha)),(1.0/(ETA_C+1.0)));
                     }
                     c2 = 0.5*((y1+y2)+betaq*(y2-y1));

                     if (c1<yl) {
                         c1=yl;
                    }
                     else if (c1>yu) {
                         c1=yu;
                     }
                     if (c2<yl){
                         c2=yl;
                     }
                     else if (c2>yu) {
                        c2=yu;
                    }

                     if (rand->nextDouble(0, 1) <= 0.5) {
                         child1->vars[i] = c2;
                         child2->vars[i] = c1;
                     }
                     else {
                         child1->vars[i] = c1;
                         child2->vars[i] = c2;
                     }
                      // }
                 }
                 else {
                     child1->vars[i] = parent1->vars[i];
                     child2->vars[i] = parent2->vars[i];
                 }
            }
            else {
                child1->vars[i] = parent1->vars[i];
                child2->vars[i] = parent2->vars[i];
            }
            
        
        }
    }
    else {
        for (int j=0; j < NUM_PARS; j++) {
                child1->vars[j] = parent1->vars[j];
                child2->vars[j] = parent2->vars[j];
        }
    }
}


void Population::simulate() {
        /* Parallel implementation */
    boost::thread *w1;
    boost::thread_group g1, g2, g3, g4, g5, g6, g7, g8;
     /* 128 threads each running 2 simulations each = 256 simulations*/
    for (int i = 0, j = 0; i < NUM_THREADS; i++, j+=NUM_SIMS) {
    
      Worker w(i, &lock, j, j+NUM_SIMS, this);
 
      w1 =  new boost::thread(w);

      if (i >=0 && i <16) {
     	g1.add_thread(w1);
      }
      else if (i >=16 && i <32) {
	g2.add_thread(w1);
      }
      else if (i >=32 && i <48) {
	g3.add_thread(w1);
      }
      else if (i >=48 && i <64){
	g4.add_thread(w1);
      }
      else if (i >=64 && i < 80) {
	g5.add_thread(w1);
      }
      else if (i >=80 && i < 96) {
	g6.add_thread(w1);
      }
      else if (i >=96 && i < 112) {
	g7.add_thread(w1);
      }
      else {
	g8.add_thread(w1);
      }
    }
    g1.join_all();
    g2.join_all();
    g3.join_all(); 
    g4.join_all();
    g5.join_all(); 
    g6.join_all();
    g7.join_all(); 
    g8.join_all();
}

void Population::merge(Population *pop1, Population *pop2) {
        /* copy parent population into positions 1-popsize */
        /* remember mixed pop is twice the size as parent and child pop */
    for (int i = 0; i < pop1->_popsize; i++) {
    //     /* copy individuals into new mixed population */
        ind[i]->copy(pop1->ind[i]);
             //ind[i] = pop1->ind[i];
        
    }     
        /*copy child population into positions popsize-2*popsize */
    for (int j = 0, k = pop2->_popsize; k < _popsize; j++, k++) {
        ind[k]->copy( pop2->ind[j]);
               // ind[k] = pop2->ind[j];
    }
}

void Population::mutate_pop() {
    for (int i=0; i<_popsize; i++) {
            /* need to pass in lower and upper limit to polynomial mutation function */
        ind[i]->mutate(rand, min_var, max_var);
    }
}


void Population::generateNewParent(Population *parent) {
    int frontCounter = 1;
    int j = 0;

    int currentSize = 0;
    
    while(frontCounter < numFronts) {

        vector<Individual *> front;
        front.clear();
        while(j < parent->_popsize) {
            if(ind[j]->rank == frontCounter) {
                front.push_back(ind[j]);
                j++;
            }
            else {
                frontCounter++;
                break;
            }
        }
        if (currentSize + (int) front.size() < parent->_popsize) {
                //cout << "Front size is " << currentSize << "\n";
            for (int p = currentSize, q =0 ; q < (int) front.size(); q++, p++) {
                parent->ind[p]->copy(front[q]);
            }
            
            currentSize += front.size(); 
        }
        else {
            sort(front.begin(), front.end(), compareDistance);
    
            for (int m =0, n = currentSize; n < parent->_popsize; m++, n++) {
                parent->ind[n]->copy(front[m]);
            }
            break;
        }
        
    }
    
}

void Population::readLimits(ifstream & ifs) {
  string line;
     double *buff = new double[2];
     if (ifs.is_open()) {  
           int i = 0;
           while (i < 8) {  
               getline( ifs, line);
               stringstream sstr( line );
               double value;
               int j = 0;
             
               while( sstr >> value) {
                   buff[j] =  value;
                   j++;          
               }
	          cout << buff[0] << "\t" << buff[1] << endl;
                  //cout << i << " " << j << endl;
              
              min_var[i] = buff[0];
              max_var[i] = buff[1];
              i++;
          }
          delete [] buff;
      }
     
    //  else {
    //      cout << "Unable to open limits file" << endl;
    //  }
}


int Population::getPopSize() {
    return _popsize;
}
