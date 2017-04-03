#ifndef __WORKER_H
#define __WORKER_H

#include "population.h"

/* even though Worker is treated as an inner nested class it does not
 * belong ot any class therefore we pass a reference to the instance
 * of the outer population
 */
class Population;

/*Runnable interface should be used if you are only planning to
 * override the run() method and no other Thread methods */
class Worker {
      private:
  //Poco::RWLock * _lock;
	boost::mutex *_mutex;
        int _id;
        int _begin;
        int _end;
        Population *_pop;
      public:
	Worker(int id, boost::mutex *mutex, int begin, int end, Population *pop);
        //Worker(int id, Poco::RWLock *lock, int begin, int end, Population *pop);
            /* When an object implementing interface Runnable is used
             * to create a thread, starting the thread causes the
             * object's run method to be called in that separately
             * executing thread. */
        //void run();
	~Worker();
	void operator()();
};

#endif
