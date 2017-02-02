#ifndef _SPENDING_RULE_H_
#define _SPENDING_RULE_H_

#include <string>
#include <functional>


class GeometricRule: public std::unary_function<double,double>
{
  const double mPsi;

  public:

  GeometricRule (double psi): mPsi(psi) { }   

  std::string identifier()    const; 
  double operator()(double w) const; 
};
  

class UniversalRule: public std::unary_function<double,double>
{
  
  public:

  UniversalRule ()  { }
  
  std::string identifier()    const;
  double operator()(double w) const;
};

#endif
