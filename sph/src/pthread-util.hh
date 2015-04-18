#ifndef sph_pthread_util_hh
#define sph_pthread_util_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

/// base class for template class Thread, defines the interfaces
/// and also implements as far as possible
struct ThreadBase {
  pthread_t pthread;
  pthread_attr_t threadattr;

  virtual ~ThreadBase() {}

  virtual void run() { }

  virtual int cancel() {
    return pthread_cancel(pthread);
  }

  virtual int join(void **returnValuePtr = 0) {
    return pthread_join(pthread, returnValuePtr);
  }

};


/// Template class Thread takes in the constructor a function pointer
/// or object of type Fun and an Argument of type Arg. The member run
/// will create a pthread which runs the function Fun with Arg as
/// argument.

/// \tparam Fun function pointer or object
/// \tparam Arg argument to pass to the function when running the thread

template<class Fun, class Arg>
struct Thread : public ThreadBase {
  Fun fun;
  Arg arg;

  Thread(Fun const &f,  Arg const &a) : 
    fun(f), arg(a) 
  {
    pthread_attr_init(&threadattr);
    int detachState = 0;
    pthread_attr_getdetachstate(&threadattr, &detachState);
    assert(detachState == PTHREAD_CREATE_JOINABLE);
  }

  ~Thread() {
    pthread_attr_destroy(&threadattr);
  }
  
  static void *runThread(void *arg) {
    Thread *thread = (Thread*) arg;
    thread->fun(thread->arg);
    return 0;
  }

  void run() {
    int const cr = pthread_create(&pthread, &threadattr, Thread::runThread, this);
    if (cr != 0) {
      std::cerr << "error: pthread_create has failed with code " << cr << "\n";
      exit(-1);
    }
  }
};

template<class Fun, class Arg>
inline Thread<Fun, Arg> *makeThread(Fun const &f, Arg const &a) {
  return new Thread<Fun,Arg>(f, a);
}

/// Utility object that wraps a pthread mutex in a C++ class.
class Mutex {
  pthread_mutex_t mutex;
  pthread_mutexattr_t mutexattr;
  bool m_locked;
public:
  Mutex(bool recursive = false) : m_locked() {
    pthread_mutexattr_init(&mutexattr);
    if (recursive) {
      pthread_mutexattr_settype(&mutexattr, PTHREAD_MUTEX_RECURSIVE);
    } else {
      pthread_mutexattr_settype(&mutexattr, PTHREAD_MUTEX_ERRORCHECK);
    }
    pthread_mutex_init(&mutex, &mutexattr);
  }
  ~Mutex() {
    pthread_mutex_destroy(&mutex);
    pthread_mutexattr_destroy(&mutexattr);
  }

  bool locked() const { return m_locked; }
  void lock() {
    pthread_mutex_lock(&mutex);
    m_locked = 1;
  }
  void unlock() {
    pthread_mutex_unlock(&mutex);
    m_locked = 0;
  }
  friend class Condition;
};

/// Utility object that wraps a pthread condition in a C++ class.
class Condition {
//   Mutex mutex;
  pthread_cond_t cond;
  pthread_condattr_t condattr;
public:
  Condition() {
    pthread_condattr_init(&condattr);
    pthread_cond_init(&cond, &condattr);
  }
  ~Condition() {
    pthread_cond_destroy(&cond);
    pthread_condattr_destroy(&condattr);
  }

  int wait(Mutex &mutex) {
    assert(mutex.locked());
    int const wr = pthread_cond_wait(&cond, &mutex.mutex);
    mutex.m_locked = false;
    return wr;
  }
  int broadcast() {
    return pthread_cond_broadcast(&cond);
  }
  int signal() {
    return pthread_cond_signal(&cond);
  }
};

#endif
