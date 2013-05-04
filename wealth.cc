#include "wealth.h"
#include "line_search.Template.h"

#include <utility>
#include <assert.h>
#include <math.h>
#include <sstream>

static const std::string messageTag ("WLTH: ");

//     DualWealthArray     DualWealthArray     DualWealthArray     DualWealthArray     DualWealthArray

void
print_pair_to_os (std::ostream& os, std::pair<double,double> const& p)
{
  os << "<" << p.first << "," << p.second << ">";
}
  

void
DualWealthArray::print_to (std::ostream& os) const
{
  os << "Dual Wealth Array '" << mName << "'[dim=" << mWealthBid.size() << "] with wealth W[0]=";
  print_pair_to_os (os, mWealthBid[0]);
  os << " and at zero index W[" << mZeroIndex << "]=";
  print_pair_to_os (os, mWealthBid[mZeroIndex]);
  os << std::endl;
}
  
void
DualWealthArray::write_to(std::ostream& os, bool asLines) const
{
  print_to(os);
  print_arrays(os, asLines);
}

void
DualWealthArray::print_arrays (std::ostream& os, bool asLines) const
{
  if (asLines)
  { os << "Wealths:  ";
    for(int i=0; i<(int) mWealthBid.size(); ++i)          // wealths
      os << mWealthBid[i].first << " ";
    os << std::endl << "  Bids: " ;
    for(int i=0; i<(int) mWealthBid.size(); ++i)          // bids
      os << mWealthBid[i].second << " ";
    os << std::endl << "Reject pos: ";
    for(int i=0; i<(int) mRejectPositions.size(); ++i)    // reject positions
      os << mRejectPositions[i].first << " ";
    os << std::endl << "Reject wt : ";
    for(int i=0; i<(int) mRejectPositions.size(); ++i)
      os << mRejectPositions[i].second << " ";
    os << std::endl << "   Bid pos: ";
    for(int i=0; i<(int) mBidPositions.size(); ++i)       // bid positions
      os << mBidPositions[i].first << " ";
    os << std::endl << "   Bid wt : ";
    for(int i=0; i<(int) mBidPositions.size(); ++i)
      os << mBidPositions[i].second << " ";
    os << std::endl;
  } else
  { os << "      <Wth,Bid>    <Rej i, p>   <Bid i, p> \n";
    for (int i=0; i<(int)mRejectPositions.size(); ++i)
    { os << " i=" << i << "   ";
      print_pair_to_os(os, mWealthBid[i]);
      os << "    ";
      print_pair_to_os(os, mRejectPositions[i]);
      os << "    ";
      print_pair_to_os(os, mBidPositions[i]);    
      os << std::endl;
    }
  }
}


std::pair<int,double>
DualWealthArray::find_wealth_position(int k, int k1, double w) const
{
  int k0 (k);
  while (k0+1 < k1)
  { int kk = floor((k0+k1)/2);
    //    std::cout << "DEBUG: init reject array   W[" << k1 << "]=" << wealth(k1)
    //	      << "  <= W[" << kk << "]=" << wealth(kk)
    //        << "  <= W[" << k0 << "]=" << wealth(k0) << std::endl;
    if (wealth(kk) < w) { k1 = kk; }
    else                { k0 = kk; }
  }
  if (k0<k1)  // inside range
  { double p ( (wealth(k0) - w) / (wealth(k0) - wealth(k1)) );
    if (p < 0) std::clog << messageTag << "*** Error ***  At k=" << k << " for w=" << w << "; share p=" << p << " to position k0=" << k0 << std::endl;
    return (std::make_pair(k0,1-p));
  }
  else return (std::make_pair(k0,1));
}
 

void
DualWealthArray::initialize_reject_array ()         
{
  double maxWealth (mWealthBid[0].first);
  mRejectPositions.push_back( std::make_pair(0,1) );                 // if reject from state 0, stay in state 0
  for (int k=1; k<(int)mWealthBid.size(); ++k)
  { double changeInWealth = mOmega - bid(k);
    double wealthAfterReject = wealth(k) + changeInWealth;
    if (maxWealth < wealthAfterReject)                               // truncate to top wealth
      mRejectPositions.push_back( std::make_pair(0,1) );
    else
    { int k0 = (changeInWealth > 0) ? 0 : k;                         // W[k1] <= new wealth < W[k0]
      int k1 = (changeInWealth > 0) ? k : mWealthBid.size()-1;       //      k0 <= k1
      mRejectPositions.push_back( find_wealth_position(k0, k1, wealthAfterReject) );
    }
  }
}

void
DualWealthArray::initialize_bid_array ()         
{
  for (int k=0; k<(int)mWealthBid.size()-1; ++k)            // not the last one
  { double wealthAfterBid = wealth(k) - bid(k);
    mBidPositions.push_back( find_wealth_position(k, mWealthBid.size()-1, wealthAfterBid) );
  }
  // fix for last bid position
  mBidPositions.push_back( std::make_pair((int)mBidPositions.size(),1) );
}

// set spacing for grid
double
DualWealthArray::grid_delta (double wealth) const
{
  if      (2.0  < wealth) return 0.2;
  else if (1.0  < wealth) return 0.1;
  else if (0.5  < wealth) return 0.05;
  else if (0.1  < wealth) return 0.01;
  else if (0.06 < wealth) return 0.005;
  else if (0.02 < wealth) return 0.002;
  else return 0.001;
}
