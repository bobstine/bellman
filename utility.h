#ifndef _UTILITY_H_
#define _UTILITY_H_

#include "wealth.h"

#include <Eigen/Core>

#include <utility>         // pair
#include <functional>
#include <math.h>
#include <assert.h>
#include <iostream>      // debug


typedef Eigen::MatrixXf Matrix;

typedef std::pair<int,double> WIndex;

////////////////////////////////////  Utility functions  /////////////////////////////////////////

double   reject_prob(double mu, double level);

double   reject_value(   int i         , WIndex const& kp , Matrix const& value, bool show = false);
double   reject_value(WIndex const& kp ,       int j      , Matrix const& value, bool show = false);
double   reject_value(WIndex const& kp1, WIndex const& kp2, Matrix const& value, bool show = false);

double   z_alpha       (double a);

double   risk          (double mu, double alpha); 

double   optimal_alpha (double mu, double omega);

////   Vector utility trades off between two possible values

class VectorUtility: public std::unary_function<double,double>
{
 protected:
  const double mAngle, mSin, mCos;
  const double mAlpha;
  double mBeta;
  double mRejectValue, mNoRejectValue;
  
 public:

 VectorUtility(double angle, double alphaLevel)
   : mAngle(angle), mSin(sin(angle * 3.1415926536/180)), mCos(cos(angle * 3.1415926536/180)),
    mAlpha(alphaLevel),  mBeta(0.0), mRejectValue(0.0), mNoRejectValue(0.0) { }
  
  double alpha      () const { return mAlpha; }
  double beta       () const { return mBeta;  }
  double angle      () const { return mAngle; }
  
  void set_constants (double beta, double rejectValue, double noRejectValue);

  double r_mu_alpha (double mu) const;
  double r_mu_beta  (double mu) const;
  
  std::pair<double,double> reject_probabilities (double mu) const;    // prob rejecting for alpha and beta
  
  virtual
    double operator()(double mu) const  { std::cout << "UTIL:  Call to operator of base class." << std::endl; return 0*mu; }

  virtual
    double bidder_utility (double mu, double rejectValue, double noRejectValue) const = 0;
  
  virtual
    double oracle_utility (double mu, double rejectValue, double noRejectValue) const = 0;
  
}; 



////  RejectVectorUtility     Rejects     Rejects     Rejects     Rejects     Rejects     Rejects     

class RejectVectorUtility: public VectorUtility

{
 public:

 RejectVectorUtility(double angle, double alpha)
   : VectorUtility(angle, alpha) { }

  double operator()(double mu) const;

  double bidder_utility (double mu, double rejectValue, double noRejectValue) const;
  double oracle_utility (double mu, double rejectValue, double noRejectValue) const;
  
}; 


////  RiskVectorUtility     Risk     Risk     Risk     Risk     Risk     Risk     Risk     Risk     Risk

class RiskVectorUtility: public VectorUtility
{
  double (*mOracleRisk)(double);
  
 public:
  
 RiskVectorUtility(double angle, double alpha)
   : VectorUtility(angle, alpha) { set_oracle_risk(); print_type(); }
  
  double operator()(double mu) const;
  
  double bidder_utility (double mu, double rejectValue, double noRejectValue) const;
  double oracle_utility (double mu, double rejectValue, double noRejectValue) const;

 private:
  void set_oracle_risk();
  void print_type() const;
  
};


//     Criteria     Criteria     Criteria     Criteria     Criteria     Criteria     Criteria     Criteria     Criteria

class AngleCriterion:  public std::binary_function<double,double,double>
{
 private:
  const double mAngle, mSin, mCos;

 public:
  AngleCriterion(double angle)
    : mAngle(angle), mSin(sin(angle * 3.1415926536/180)), mCos(cos(angle * 3.1415926536/180)) { }

  std::string identifier() const     { return std::to_string(mAngle); }

  double operator()(double x, double y) const { return mCos*x + mSin*y; }
};


class RiskInflationCriterion:  public std::binary_function<double,double,double>
{
 private:
  const double mB1;

 public:
  RiskInflationCriterion(double b1)
    : mB1(b1) { }

  std::string identifier() const     { return std::to_string((int)round(mB1)); }

  double operator()(double x, double y) const { return x - mB1*y; }
};


//  Matrix     Matrix     Matrix     Matrix     Matrix     Matrix     Matrix     Matrix     Matrix     Matrix     Matrix

class MatrixUtility: public std::unary_function<double,double>
{
 protected:
  double mAlpha, mBeta;
  double mV00, mV01, mV10, mV11;      // 0 for not reject, 1 for reject

 public:
  
 MatrixUtility()
   : mV00(0.0), mV01(0.0), mV10(0.0), mV11(0.0) {}
  
  void set_constants (double alpha, double beta, double v00, double v01, double v10, double v11);
  
  double r_mu       (double mu, double p) const;
  double r_mu_alpha (double mu)           const;
  double r_mu_beta  (double mu)           const;
  
  std::pair<double,double> reject_probabilities (double mu) const;    // prob rejecting for alpha and beta
}; 



////  RejectMatrixUtility     Rejects     Rejects     Rejects     Rejects     Rejects     Rejects     

template <class C>
class RejectMatrixUtility: public MatrixUtility
{
 private:
  C mCriterion;

 public:
  
 RejectMatrixUtility(C criterion)
   : MatrixUtility(), mCriterion(criterion) { }

  std::string identifier() const                { return mCriterion.identifier(); }

  double operator()(double mu) const;

  double row_utility (double mu, double v00, double v01, double v10, double v11) const;
  double col_utility (double mu, double v00, double v01, double v10, double v11) const;  

  void   write_details_to_stream(double mu, std::ostream& os) const;
}; 


////  Risk     Risk     Risk     Risk     Risk     Risk     Risk     Risk     Risk     Risk

template <class C>
class RiskMatrixUtility: public MatrixUtility
{
 private:
  C mCriterion;

 public:

 RiskMatrixUtility(C criterion)
   : MatrixUtility(), mCriterion(criterion) { }
  
  std::string identifier() const                { return mCriterion.identifier(); }
  
  double operator()(double mu) const;
  
  double negative_risk(double mu, double alpha) const;

  double row_utility (double mu, double v00, double v01, double v10, double v11) const;
  double col_utility (double mu, double v00, double v01, double v10, double v11) const;

  void   write_details_to_stream(double mu, std::ostream& os) const;
  
}; 



#endif
