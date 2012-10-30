#ifndef _DISTRIBUTION_H_
#define _DISTRIBUTION_H_

#include <string>
#include <functional>

class Distribution: public std::unary_function<int,double>
{
  public:
  virtual
    std::string identifier()    const = 0;
  virtual
    double      operator()(int) const = 0;
};



/**********************************************************************************

   Probability distributions that control spending

      Zero origin: All operator() defined so that \sum_{k=0} p(k) = 1 

 **********************************************************************************/

class GeometricDist: public Distribution
{
  const double mPsi;
  const double mOneMinusPsi;

  public:

  GeometricDist (double psi): mPsi(psi), mOneMinusPsi(1-psi) { }   

  std::string identifier() const; 
  double operator()(int k) const;    // percent of wealth spent at step k; adds to 1
};
  
class UniformDist: public Distribution
{
  const int     mLimit;
  const double  mP;     // 1/(number of tests)
  
 public:

  UniformDist (double n): mLimit(n), mP(1.0/n) {}
  
  std::string identifier() const;
  double operator()(int k) const;
};

class UniversalDist: public Distribution
{
  const int mStart;  // starting index
  
  public:

  UniversalDist (int start): mStart (start) { }
  
  std::string identifier() const;
  double operator()(int k) const;
  };


//  scaled     scaled     scaled     scaled     scaled     scaled     scaled

/*
  A scaled universal distribution starts with the function

      g(k) = c/(k log^2(k+1))   k = 1, 2, ...

  The scaling constant c is a bit arbitrary though one can make an argument for
  k about 4 based on a recurrance (Deans notes).  g determines the bid amounts
  for the first bid (the maximum bid, k=1), second and so forth on down.  We set
  the wealth by inserting the sum of g(k) into the wealth function at W[0] and
  decrementing that sum by g(k).  The find indexing function requires adjustment
  as well since some of the bid amounts are larger than the payout.
*/

class UniversalBidder: public std::unary_function<double,double>
{
  double const mScale;                        // k in dean's notes
  
 public:
  
  UniversalBidder (double scale)    : mScale(scale)  { }
  
  std::string identifier()                const;
  double      operator()(double round)    const;
  
  double      total_wealth()              const;
};
  


class ScaledUniversalDist: public Distribution
{
  static const double mSumOfRecipLog;
  const double mScale;                        // k in dean's notes
  
 public:
  
  ScaledUniversalDist (double scale)    : mScale(scale)  { }
  
  std::string identifier()           const;
  double operator()(int k)           const;
  
  int w0_index(double initialWealth) const;   // int such that scaled tail sum matches initial wealth
  double   max_wealth()              const  { return mScale * mSumOfRecipLog; }

 private:  
  double g(int k)                    const;   // the log recip with scaling factor
};


// double uniform_to_end (int k, int left);

#endif
