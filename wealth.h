#ifndef _WEALTH_H_
#define _WEALTH_H_

#include "dynamic_array.h"

#include <iostream>
#include <functional>
#include <vector>

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

class ScaledUniversalDist: public Distribution
{
  static const double mSumOfRecipLog;

  const double mScale;          // k in dean's notes

 public:

  ScaledUniversalDist (double scale)    : mScale(scale)  { ; }
  
  std::string identifier() const;
  double operator()(int k) const;

  double max_wealth() const { return mScale * mSumOfRecipLog; }
  
  int starting_index(double initialWealth) const;   // int such that scaled tail sum matches initial wealth
  
};

// double uniform_to_end (int k, int left);


/**********************************************************************************

   Wealth tracker basically maps from integers that follow the tests
   to the wealth of the expert or bidder.  The wealth is monotone increasing
   in the index/position k.

   Normalized to have wealth omega at index mZeroIndex.

 **********************************************************************************/

class WealthArray
{
  const std::string     mName;
  const int             mPadding;    // space for wealth above omega; keep < 15 or cannot solve for multiplier.
  const int             mZeroIndex;  // position of W_0, the place used for omega; set to max number steps can take
  const int             mSize;       // number of distinct wealth values
  const double          mOmega;      // defines wealth at zeroIndex and determines how far 'up' wealth can go 
  DynamicArray<double>  mWealth;     // indices k < mZeroIndex denote wealth less than omega
  std::vector< std::pair<int,double> > mPositions;  // cache locations for new positions when increment wealth by rejection

 public:

  WealthArray ()
    : mName("empty"), mPadding(0), mZeroIndex(0), mSize(mPadding), mOmega(0), mWealth(), mPositions() { }
  
 WealthArray(double w0, double omega, int zeroIndex, ScaledUniversalDist const& f)
   : mName(f.identifier()),   mPadding(f.starting_index(w0)), mZeroIndex(zeroIndex), mSize(zeroIndex+mPadding), mOmega(omega), mWealth(), mPositions() { initialize_array_using_func(f);}
  
 WealthArray(double omega, int zeroIndex, Distribution const& pdf)
   : mName(pdf.identifier()), mPadding(2),                    mZeroIndex(zeroIndex), mSize(zeroIndex+mPadding), mOmega(omega), mWealth(), mPositions() { initialize_array_using_pdf(pdf);}

 WealthArray(double omega, int zeroIndex, double psi) // use for geometric for numerical stability
   : mName(geom_name(psi)),   mPadding(2),                    mZeroIndex(zeroIndex), mSize(zeroIndex+mPadding), mOmega(omega), mWealth(), mPositions() { initialize_geometric_array(psi);}


  std::string name()               const { return mName; }
  int    size ()                   const { return mSize; }
  int    zero_index ()             const { return mZeroIndex ; }
  double omega ()                  const { return mOmega; }
  
  double bid(int k)                const { return mWealth[k+1]-mWealth[k]; }    // recursion goes from small wealth to larger (note inlined in init array too)
  double wealth(int k)             const { return mWealth[k]; }
  double operator[](int k)         const { return mWealth[k]; }
  
  std::pair<int, double> wealth_position (int k) const { return mPositions[k]; }
  
  void print_to (std::ostream& os) const { os << "Wealth array " << mName << "  " << mWealth; }
  
 private:
  std::string geom_name(double p) const;
  void init_positions ();
  void initialize_array_using_pdf (Distribution const& p);
  void initialize_array_using_func(ScaledUniversalDist const& p);
  void initialize_geometric_array (double psi);
  void fill_array_top();
  std::pair<int, double> find_wealth_position (int k, double increaseInWealth) const;

};

inline
std::ostream&
operator<< (std::ostream& os, WealthArray const& wa)
{
  wa.print_to(os);
  return os;
}



#endif
