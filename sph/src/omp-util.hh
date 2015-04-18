#ifndef sph_omp_util_hh
#define sph_omp_util_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

///////////////////////////////////////////////////////////////////////////////////////
/// Utility class that obtains the team-id of an OMP-Thread. The object can either be
/// printed with << or it can provide the team id as a number.			     
///////////////////////////////////////////////////////////////////////////////////////
struct OMPThread {
  /// Hold the thread's team id.
  size_t const myid;
  /// Hold the thread's team size.
  size_t const teamsize;

  /// Default-Constructor. 
  /// This initializes myid, teamsize by calls to omp_get_thread_num()
  /// omp_get_num_threads()
  OMPThread() : 
    myid(omp_get_thread_num()), 
    teamsize(omp_get_num_threads()) 
  {}

  /// Provide the team id as a number.
  /// \returns myid
  operator size_t() const { return myid; }

  /// Print the object to an std::ostream.
  /// \param aus an ostream
  void print(std::ostream &aus) const {
    aus << "Thread " << myid << "/" << teamsize;
  }
};
/// Overload of operator << to print an OMPThread to a stream.
/// \param aus the stream
/// \param p an OMPThread object
/// \return the stream
inline std::ostream &operator <<(std::ostream &aus, OMPThread const &p) {
  p.print(aus);
  return aus;
}


///////////////////////////////////////////////////////////////////////////////////////
/// Wrapper object that provides an OMP lock or mutex.			     
/// This lock does not allow for nested locking. If you need that see OMPNestLock.
///////////////////////////////////////////////////////////////////////////////////////
class OMPLock {
  /// Holds the opaque lock data.
  omp_lock_t theLock;

  /// Holds the optional name of the lock.
  std::string const name;

  /// Holds the team id of the current owner of the lock.
  size_t curOwner;

  /// Remembers whether the Lock is locked.
  bool isLocked;

public:
  /// Default constructor. 
  /// Initializes theLock by calling omp_init_lock().
  OMPLock() : curOwner(), isLocked() {
    omp_init_lock(&theLock);
  }

  /// This constructor allows to set a name for the lock. 
  /// Initializes theLock by calling omp_init_lock().
  /// \param name the lock's name
  OMPLock(std::string const &name) : 
    name(name), 
    curOwner(), 
    isLocked() {
    omp_init_lock(&theLock);
  }

  /// The destructor destroys the OMP lock by calling omp_destroy_lock() on theLock.
  ~OMPLock() {
    omp_destroy_lock(&theLock);
  }

  /// Locks the lock via omp_set_lock().
  void lock() {
    omp_set_lock(&theLock);
    isLocked = 1;
    curOwner = omp_get_thread_num();
  }

  /// Unlocks the lock via omp_unset_lock().
  void unlock() {
    isLocked = 0;
    curOwner = 0;
    omp_unset_lock(&theLock);
  }

  /// Test if lock is obtainable by calling omp_test_lock().
  /// \return result of call to omp_test_lock()
  int testlock() {
    int const r = omp_test_lock(&theLock);
    if (r) {
      isLocked = 1;
      curOwner = omp_get_thread_num();
    }
    return r;
  }

  /// Return true if the lock is locked currently.
  /// \return isLocked.
  bool locked() const {
    return isLocked;
  }
  /// Return the team id of the current owner of the lock.
  /// \return curOwner.
  size_t owner() const {
    return curOwner;
  }

  /// Print the object to an std::ostream.
  /// \param aus an ostream
  void print(std::ostream &aus) const {
    aus << "Lock " << name << "(" << locked() << ")";
  }
};
/// Overload of operator << to print an OMPLock to a stream.
/// \param aus the stream
/// \param p an OMPThread object
/// \return the stream
inline std::ostream &operator <<(std::ostream &aus, OMPLock const &p) {
  p.print(aus);
  return aus;
}

///////////////////////////////////////////////////////////////////////////////////////
/// Wrapper object that provides a nested OMP lock.			     
/// If you don't need a nested lock, see OMPLock.
///////////////////////////////////////////////////////////////////////////////////////
class OMPNestLock {
  /// Holds the opaque lock data.
  omp_nest_lock_t theLock;

  /// Holds the optional name of the lock.
  std::string const name;

  /// Holds the team id of the current owner of the lock.
  size_t curOwner;

  /// Remembers whether the Lock is locked.
  bool isLocked;

public:
  /// Default constructor. 
  /// Initializes theLock by calling omp_init_nest_lock().
  OMPNestLock() : curOwner(), isLocked() {
    omp_init_nest_lock(&theLock);
  }
  /// This constructor allows to set a name for the lock. 
  /// Initializes theLock by calling omp_init_nest_lock().
  /// \param name the lock's name
  OMPNestLock(std::string const &name) : name(name), curOwner(), isLocked() {
    omp_init_nest_lock(&theLock);
  }

  /// The destructor destroys the OMP lock by calling omp_destroy_nest_lock() on theLock.
  ~OMPNestLock() {
    omp_destroy_nest_lock(&theLock);
  }

  /// Locks the lock via omp_set_nest_lock().
  void lock() {
    omp_set_nest_lock(&theLock);
    isLocked = 1;
    curOwner = omp_get_thread_num();
  }

  /// Unlocks the lock via omp_unset_nest_lock().
  void unlock() {
    isLocked = 0;
    curOwner = 0;
    omp_unset_nest_lock(&theLock);
  }

  /// Test if lock is obtainable by calling omp_test_lock().
  /// \return result of call to omp_test_nest_lock()
  int testlock() {
    int const r = omp_test_nest_lock(&theLock);
    if (r) {
      isLocked = 1;
      curOwner = omp_get_thread_num();
    }
    return r;
  }

  /// Return true if the lock is locked currently.
  /// \return isLocked.
  bool locked() const {
    return isLocked;
  }

  /// Return the team id of the current owner of the lock.
  /// \return curOwner.
  size_t owner() const {
    return curOwner;
  }

  /// Print the object to an std::ostream.
  /// \param aus an ostream
  void print(std::ostream &aus) const {
    aus << "NestLock " << name << "(" << locked() << ")";
  }
};
/// Overload of operator << to print an OMPNestLock to a stream.
/// \param aus the stream
/// \param p an OMPNestLock object
/// \return the stream
inline std::ostream &operator <<(std::ostream &aus, OMPNestLock const &p) {
  p.print(aus);
  return aus;
}


///////////////////////////////////////////////////////////////////////////////////////
/// A threadsafe ostream for a fixed maximum number of threads.  The
/// class prints all data to a std::ostringstream. If the data in that
/// string stream contains a '\\n' character, the data is printed to
/// the intented ostream with a prefix describing the thread id. The
/// template parameter N determines the maximum number of threads
/// (i.e. the maximum thread team id) allowed.
///////////////////////////////////////////////////////////////////////////////////////
template<int N=8>
class TSOStream {
  
  /// The ostream where the data shalle be written to.
  std::ostream &aus;

  /// An OMPLock controls access to the stream aus.
  OMPLock lock;

  /// N std::ostringstreams serve as buffers for each thread's data.
  std::ostringstream buffer[N];

public:

  /// Constructor.
  /// \param aus the stream to control.
  TSOStream(std::ostream &aus) : 
    aus(aus) 
  { }

  /// Write some data to the stream, in a thread-safe way.
  /// \param v the data to write
  /// \returns reference to myself
  template<class T>
  TSOStream &operator <<(T const &v) {
    std::ostringstream &sbuffer = buffer[OMPThread()];
    sbuffer << v;
    if (sbuffer.str().find('\n') != std::string::npos) {
      lock.lock();
      aus << OMPThread() << ": " << sbuffer.str();
      lock.unlock();
      sbuffer.str("");
    }
    return *this;
  }
};

#endif
