#ifndef integr_integr_base_hh
#define integr_integr_base_hh
/*
  This file is part of Sphere.
  Copyright © 2008,2009,2010 Johannes Willkomm
  See the file sph/src/sph.hh for copying conditions.  
*/

template<class Fn>
struct IntegratorBase {

  typedef Fn ForceFunctor;
  typedef typename Fn::value_type value_type;
  typedef typename Fn::first_argument_type  first_argument_type;
  typedef typename Fn::second_argument_type second_argument_type;
  typedef typename Fn::result_type result_type;
  typedef typename Fn::result_type::value_type element_type;

  IntegratorBase() { }
  virtual ~IntegratorBase() { }

  virtual result_type operator()(first_argument_type const dt, second_argument_type const &state) = 0;
  virtual void update(second_argument_type &state, result_type const &update) = 0;
  virtual std::string name() const = 0;
  virtual std::string description() const = 0;
  virtual void printXML(std::ostream &str) const = 0;

};

template<class A>
inline std::ostream &operator <<(std::ostream &str, IntegratorBase<A> const &pk) {
  pk.printXML(str);
  return str;
}


#endif
