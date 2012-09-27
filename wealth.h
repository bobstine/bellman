#ifndef _WEALTH_H_
#define _WEALTH_H_

#include "distribution.h"

#include <vector>
#include <iostream>

/**********************************************************************************

   Wealth array maps from integers that follow the tests to the wealth
   of the expert or bidder.  The wealth is monotone *decreasing* in the
   index/position k.

   Normalized to have wealth omega at index mZeroIndex.

 **********************************************************************************/

class WealthArray
{
  const std::string     mName;
  const int             mZeroIndex;  // position of W_0, the place used for locating initial wealth W0
  const double          mOmega;      // defines wealth at zeroIndex and determines how far increases
  std::vector<double>   mWealth;     
  std::vector< std::pair<int,double> > mPositions;  // cache new positions when increment wealth by rejection

 public:

  WealthArray ()
    : mName("empty"), mZeroIndex(0), mOmega(0), mWealth(), mPositions() { }
  
 WealthArray(double w0, double omega,                int steps, ScaledUniversalDist const& f)
   : mName(f.identifier()),   mZeroIndex(f.w0_index(w0)), mOmega(omega), mWealth(mZeroIndex+steps+1), mPositions() { initialize_array_using_func(f);}
  
 WealthArray(           double omega, int zeroIndex, int steps, Distribution const& pdf)
   : mName(pdf.identifier()), mZeroIndex(zeroIndex)     , mOmega(omega), mWealth(zeroIndex+steps+1), mPositions() { initialize_array_using_pdf(pdf);}

 WealthArray(           double omega, int zeroIndex, int steps, double psi) // use for geometric for numerical stability
   : mName(geom_name(psi)),   mZeroIndex(zeroIndex)     , mOmega(omega), mWealth(zeroIndex+steps+1), mPositions() { initialize_geometric_array(psi);}


  std::string name()               const { return mName; }
  int    size ()                   const { return (int) mWealth.size(); }
  int    zero_index ()             const { return mZeroIndex ; }
  double omega ()                  const { return mOmega; }
  
  double bid(int k)                const { return mWealth[k]-mWealth.at(k+1); }
  double wealth(int k)             const { return mWealth.at(k); }
  double operator[](int k)         const { return mWealth.at(k); }
  
  std::pair<int, double> wealth_position (int k) const { return mPositions[k]; }  // access to cached array
  
  void print_to (std::ostream& os) const;

  
 private:
  std::string geom_name(double p) const;
  void init_check() const;
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
