#ifndef _SPENDING_RULE_H_
#define _SPENDING_RULE_H_

#include <string>
#include <functional>

class SpendingRule: public std::unary_function<double,double>
{
  public:
  virtual
    std::string identifier()         const = 0;
  virtual
    double      operator()(double w) const = 0;           // returns alpha(wealth)
};



class GeometricRule: public SpendingRule
{
  const double mPsi;

  public:

  GeometricDist (double psi): mPsi(psi) { }   

  std::string identifier()    const; 
  double operator()(double w) const; 
};
  

class UniversalRule: public SpendingRule
{
  
  public:

  UniversalRule ()  { }
  
  std::string identifier()    const;
  double operator()(double w) const;
};

#endif
