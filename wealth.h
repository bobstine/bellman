#ifndef _WEALTH_H_
#define _WEALTH_H_

#include "dynamic_array.h"
#include "distribution.h"

#include <iostream>
#include <vector>

/**********************************************************************************

   Wealth array maps from integers that follow the tests to the wealth
   of the expert or bidder.  The wealth is monotone *decreasing* in the
   index/position k.

   Normalized to have wealth omega at index mZeroIndex.

 **********************************************************************************/

class WealthArray
{
  const std::string     mName;
  const int             mPadding;    // space for wealth above omega; keep < 15 or cannot solve for multiplier.
  const int             mZeroIndex;  // position of W_0, the place used for omega; set to max number steps can take
  const int             mSize;       // number of distinct wealth values
  const double          mOmega;      // defines wealth at zeroIndex and determines how far the wealth can increase
  DynamicArray<double>  mWealth;     // indices k < mZeroIndex denote wealth less than omega
  std::vector< std::pair<int,double> > mPositions;  // cache locations for new positions when increment wealth by rejection

 public:

  WealthArray ()
    : mName("empty"), mPadding(0), mZeroIndex(0), mSize(mPadding), mOmega(0), mWealth(), mPositions() { }
  
 WealthArray(double w0, double omega, int zeroIndex, ScaledUniversalDist const& f)
   : mName(f.identifier()),   mPadding(f.w0_index(w0)), mZeroIndex(zeroIndex), mSize(zeroIndex+mPadding), mOmega(omega), mWealth(), mPositions() { initialize_array_using_func(f);}
  
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
  
  void print_to (std::ostream& os) const;

  
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
