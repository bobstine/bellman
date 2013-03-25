#ifndef _WEALTH_H_
#define _WEALTH_H_

#include "distribution.h"

#include <vector>
#include <iostream>

/**********************************************************************************

   Wealth array maps from integers that follow the tests to the wealth
   of the expert or bidder.  The wealth is monotone *decreasing* in the
   index/position k.  Bid is difference in wealths.  This code mostly replaced
   by the DualWealthArray below.

   Normalized to have wealth omega at index mZeroIndex.

 **********************************************************************************/

class WealthArray
{
  const std::string     mName;
  const int             mZeroIndex;  // position of W_0, the place used for locating initial wealth W0
  const double          mOmega;      // defines wealth at zeroIndex; max wealth is set by zero index in fill_array_top
  std::vector<double>   mWealth;     
  std::vector< std::pair<int,double> > mPositions;  // cache new positions when increment wealth by rejection

 public:

  WealthArray ()
    : mName("empty"), mZeroIndex(0), mOmega(0), mWealth(), mPositions() { }
  
 WealthArray(double w0, double omega,                int steps, ScaledUniversalDist const& f)
   : mName(f.identifier()),   mZeroIndex(f.w0_index(w0)), mOmega(omega), mWealth(mZeroIndex+steps+1), mPositions() { initialize_array_using_func(f);}
  
 WealthArray(           double omega, int zeroIndex, int steps, Distribution const& pdf)
   : mName(pdf.identifier()), mZeroIndex(zeroIndex)     , mOmega(omega), mWealth(mZeroIndex+steps+1), mPositions() { initialize_array_using_pdf(pdf);}

 WealthArray(           double omega, int zeroIndex, int steps, double psi) // use for geometric for numerical stability
   : mName(geom_name(psi)),   mZeroIndex(zeroIndex)     , mOmega(omega), mWealth(mZeroIndex+steps+1), mPositions() { initialize_geometric_array(psi);}


  std::string name()               const { return mName; }
  int    number_of_bids()          const { return (int) mWealth.size() - 1; } // -1 to avoid bidding off end of array
  int    zero_index ()             const { return mZeroIndex ; }
  double omega ()                  const { return mOmega; }
  
  double wealth(int k)             const { return mWealth.at(k); }
  int    reject_jumps_to(int k)    const { return mPositions.at(k).first; }
  double reject_jump_share(int k)  const { return mPositions.at(k).second; }
  
  double operator[](int k)         const { return mWealth.at(k); }
  double bid(int k)                const;
  
  std::pair<int, double> wealth_position (int k) const { return mPositions[k]; }  // access to cached array
  
  void print_to (std::ostream& os) const;
  void write_to (std::ostream& os) const;   // more details
  
 private:
  std::string geom_name(double p) const;
  int  size()       const { return mWealth.size(); }
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




// Dual wealth array randomizes when reject or accept H0
// Discrete grid along wealth axis, with interpolation to find
// reject and bid fails positions.

class DualWealthArray
{
  typedef std::pair<int, double>  Position;
  typedef std::vector< Position > PositionVector;
  
  const std::string     mName;
  const double          mW0;                          // initial wealth
  double const          mOmega;                       // payout when reject
  int                   mZeroIndex;                   // position of W_0, the place used for locating initial wealth W0
  std::vector< std::pair<double,double> > mWealthBid;     
  PositionVector        mRejectPositions;             // cache new <index,prob> positions when increment wealth by rejection
  PositionVector        mBidPositions;                //       new positions when decrement wealth by bid

 public:

  DualWealthArray ()
    : mName("empty"), mW0(0), mOmega(0), mZeroIndex(0), mWealthBid(), mRejectPositions(), mBidPositions() { }

  DualWealthArray (double w0)      // constant alpha bidder with no wealth constraint; always bids w0
    : mName("fixed"), mW0(w0), mOmega(0), mZeroIndex(0), mWealthBid(1), mRejectPositions(1), mBidPositions(1)
    { mWealthBid[0] = std::make_pair(w0,w0); mRejectPositions[0] = std::make_pair(0,1); mBidPositions[0]=std::make_pair(0,1); }
  
  template <class Bidder>
    DualWealthArray(std::string name, double w0, double omega, Bidder const& f, int nRounds)
    : mName(name), mW0(w0), mOmega(omega), mZeroIndex(0), mWealthBid(), mRejectPositions(), mBidPositions()
    { initialize_wealth_bid_array(f, nRounds); initialize_reject_array(); initialize_bid_array(); }

  std::string name()               const { return mName; }
  double      initial_wealth()     const { return mW0; }
  double      omega ()             const { return mOmega; }
  int         zero_index ()        const { return mZeroIndex ; }
  int  number_wealth_positions()   const { return (int) mWealthBid.size(); }
  
  double wealth(int k)             const { return mWealthBid.at(k).first; }
  double bid (int k)               const { return mWealthBid.at(k).second; }

  Position reject_position(int k)  const { return mRejectPositions.at(k); }
  int      reject_index (int k)    const { return mRejectPositions.at(k).first; }
  double   reject_prob (int k)     const { return mRejectPositions.at(k).second; }

  Position bid_position (int k)    const { return mBidPositions.at(k); }
  int    bid_jumps_to(int k)       const { return mBidPositions.at(k).first; }
  double bid_jump_share(int k)     const { return mBidPositions.at(k).second; }
  
  std::pair<double,double> operator[](int k)  const { return mWealthBid.at(k); }
  

  void print_to (std::ostream& os) const;
  void write_to (std::ostream& os, bool asLines=false) const;   // more details
  
 private:
  template <class Bidder>
  void initialize_wealth_bid_array (Bidder const& f, int nRounds);
  void initialize_reject_array ();
  void initialize_bid_array    ();
  std::pair<int,double> find_wealth_position(int loIndex, int hiIndex, double w) const;
  double grid_delta (double wealth) const;
  void print_arrays   (std::ostream& os, bool asLines=false) const;
};




inline
std::ostream&
operator<< (std::ostream& os, DualWealthArray const& wa)
{
  wa.print_to(os);
  return os;
}



#endif
