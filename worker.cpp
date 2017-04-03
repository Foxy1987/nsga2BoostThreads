#include "worker.h"


    /* each thread will be assigned a portion of the
     * individual array to simulate */
Worker::Worker(int id, boost::mutex* mutex, int begin, int end, Population *pop) :
    _id(id), _mutex(mutex), _begin(begin), _end(end), _pop(pop) {}

Worker::~Worker()  {
}

void Worker::operator()() {
    for (int i = _begin; i < _end; i++) {
            /* Current thread will simulate child strings from begin
             * to end */
	/* When the lock is constructed it holds onto the mutex until
	   it is destructed. */ 
	//boost::mutex::scoped_lock _lock(_mutex); 
	//boost::mutex::scoped_lock scoped_lock(*_mutex);
	        Individual *ind = _pop->ind[i];
        ind->simulate(_id);
	_mutex->lock(); //acquire lock



        //cout << "Thread " << _id << " is simulating child at position " << i << "\n";
	_mutex->unlock();  //release lock
    }
    
}



// #include "worker.h"


//     /* each thread will be assigned a portion of the
//      * individual array to simulate */
   
// Worker::Worker(int id, Poco::RWLock *lock, int begin, int end, Population *pop) 
//     : _id(id), _lock(lock), _begin(begin), _end(end), _pop(pop) {}

// void Worker::run() {
//     for (int i = _begin; i < _end; i++) {
//             /* Current thread will simulate child strings from begin
//              * to end */
//         Individual *ind = _pop->ind[i];
//             /* the simulate method runs the simulator, zfgenerator and
//              * saves stats in individual members */
//         ind->simulate(_id);
                
//         _lock->writeLock(); //acquire lock
//         cout << "Thread " << _id << " is simulating child at position " << i << "\n";
//         _lock->unlock();  //release lock
//     }
// }
