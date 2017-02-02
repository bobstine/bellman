#ifndef _WEALTH_H_
#define _WEALTH_H_

#include "spending_rule.h"

#include <vector>
#include <iostream>


// Dual wealth array randomizes when reject or accept H0
// Puts a discrete grid along wealth axis, with interpolation to find reject and bid fails positions.

class DualWealthArray
{
  typedef std::pair<int, double>  Position;           // first is index, second is probability (with rest on next higher index wealth with 1-p)
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
  
  template <class Rule>
    DualWealthArray(std::string name, double maxWealth, double w0, double omega, Rule const& f, int nRounds)
    : mName(name), mW0(w0), mOmega(omega), mZeroIndex(0), mWealthBid(), mRejectPositions(), mBidPositions()
    { initialize_wealth_bid_array(f, maxWealth, nRounds); initialize_reject_array(); initialize_bid_array(); }

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
  template <class Rule>
    void initialize_wealth_bid_array (Rule const& f, double maxWealth, int nRounds);
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
