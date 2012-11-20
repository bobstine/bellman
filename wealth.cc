#include "wealth.h"
#include "line_search.Template.h"

#include <utility>
#include <assert.h>
#include <math.h>
#include <sstream>

static const std::string messageTag ("WLTH: ");


double
WealthArray::bid(int k)                const
{
  assert ((0 <= k) && (k < size()-1));   // cannot find bid past end of array
  double b = mWealth[k]-mWealth[k+1];
  if (1.0 < b)
  { std::clog << messageTag << "Wealth array bid at k=" << k << " is larger than 1 (" << b << "; reduced to 1." << std::endl;
    b = 1;
  }
  return b;
}

std::pair<int, double>
WealthArray::find_wealth_position (int k, double increase)  const // k0 is current position denoting current wealth
{
  int k1 (k);
  double target = mWealth[k1] + increase;       // 'wealth' is 'new wealth' > 'current wealth'
  int k0 (0);                            
  while (mWealth[k1] > target)                  // if bid > payoff
    ++k1;
  while (k0+1 < k1)                             // bracket target by wealth[k0] >= target >= wealth[k1]
  { int kk = floor((k0+k1)/2);
    // std::cout << "DEBUG: find_wealth_position   W[" << k0 << "]="
    //           << mWealth[k0] << " W[" << kk << "]=" << mWealth[kk] << "  W[" << k1 << "]=" << mWealth[k1] << std::endl;
    if (mWealth[kk] < target)
    { k1 = kk; }
    else
    { k0 = kk; }
  }
  if (k0<k1)  // inside range
  { double p ( (mWealth[k0] - target) / (mWealth[k0] - mWealth[k1]) );
    if (p < 0) std::cerr << messageTag << "*** Error ***  At k=" << k << " for increment " << increase << "; share p=" << p << " to position k0=" << k0 << std::endl;
    return std::make_pair(k0,1-p);
  }
  else
    return std::make_pair(k0,1);
}


void
WealthArray::init_check() const
{
  std::clog << messageTag << "Initializing wealth array with " << size() << " values for run with "<< size()-mZeroIndex
	    << " steps and wealth " << mOmega << " @ " << mZeroIndex << std::endl;
  assert((0 < mZeroIndex) && (mZeroIndex < size()));
}

void
WealthArray::initialize_array_using_pdf(Distribution const& p)
{
  init_check();
  mWealth[mZeroIndex]=mOmega;
  for(int i=mZeroIndex+1; i < size(); ++i)
  { double bid (mOmega * p(i-mZeroIndex-1));
    // std::cout << "    wealth[i="<<i<<"] = (wealth["<<i-1<<"]="<< mWealth[i-1]<<")-("<<mOmega<<")*(p["<< i-mZeroIndex-1 <<"]="<< p(i-mZeroIndex-1)<<")\n";
    mWealth[i] = mWealth[i-1] - bid;
  }
  fill_array_top();
  init_positions();
}


void
WealthArray::initialize_array_using_func(ScaledUniversalDist const& f)
{
  init_check();
  mWealth[0] = f.max_wealth();
  std::clog << messageTag << "Initial wealth W[0]=" << mWealth[0] << std::endl;
  for(int i=1; i<size(); ++i)
    mWealth[i] = mWealth[i-1]-f(i-1);
  init_positions();
}


void
WealthArray::initialize_geometric_array(double psi)
{
  init_check();
  mWealth[mZeroIndex] = mOmega;
  for(int i=mZeroIndex+1; i < size(); ++i)
  { double bid (mWealth[i-1]*psi);
    mWealth[i] = mWealth[i-1] - bid; 
  }
  fill_array_top();
  init_positions();
}


//  fill wealth above omega by incrementing omega over iZero steps
void
WealthArray::fill_array_top()
{
  double w (0.5);                   // allow to grow by this much 
  int    k (mZeroIndex) ;           // over this many steps
  // geometric sum
  double b (bid(mZeroIndex));       // incrementing initial bid
  double m (Line_Search::Bisection(0.00001,std::make_pair(1.000001,3))
	    ([&w,&k,&b](double x){ double xk(x); for(int j=1;j<k;++j) xk *= x; return x*(1.0-xk)/(1-x) - w/b;}));
  std::cout << messageTag << "Geometric to top, b=" << b << " and growth factor m= " << m << " applied to " << mWealth[mZeroIndex] <<std::endl;
  if (m < 1)
  { m = 1.0;
    std::cerr << messageTag << " *** Error ***  Wealth array cannot initialize upper wealth for inputs. Setting m = 1." << std::endl;
    std::cerr << "            w=" << w << "    k=" << k << "   b=" << b << std::endl;
  }
  for(int i=mZeroIndex-1; 0 <= i; --i)
  { b *= m;
    mWealth[i] = mWealth[i+1] + b;
  }
  // check that top bid is at least omega
  if (bid(0) < mOmega)
  { std::clog << messageTag << "Note: Moving bid at zero to omega."<< std::endl;
    mWealth[0] = mWealth[1]+mOmega;
  }
}


void
WealthArray::init_positions ()
{
  // cache indexing for new positions since the increment omega is fixed; -1 since do not have 'next' at end.
  for(int j = 0; j<size()-1; ++j)  
  { double increase = mOmega - bid(j);
    if (increase < 0)
      std::cerr << messageTag << "*Warning*  Wealth implies certain loss because bid " << bid(j) << " exceeds payoff " << mOmega  << std::endl;
    mPositions.push_back( find_wealth_position(j,increase) );
  }
  //  for (int i=0; i<size()-1; ++i) 
  //  std::cout << messageTag << "DEBUG i=" << i << "  $" <<mWealth[i] << "  " << bid(i)
  //	      << "  [" << mPositions[i].first << "," << mPositions[i].second << "]" << std::endl;
}


void
WealthArray::print_to (std::ostream& os) const
{
  os << "Wealth array " << mName << "[dim=" << mWealth.size() << "] with wealth vector beginning W[0]=" << mWealth[0]
     << " and at zero index W[" << mZeroIndex << "]=" <<  mWealth[mZeroIndex] << std::endl;
}
  
void
WealthArray::write_to(std::ostream& os) const
{
  for (int i=0; i<size(); ++i)               os << mWealth[i]           << " ";    os << std::endl;
  for (int i=0; i<number_of_bids(); ++i)     os << reject_jumps_to(i)   << " ";    os << std::endl;
  for (int i=0; i<number_of_bids(); ++i)     os << reject_jump_share(i) << " ";    os << std::endl;
}


std::string
WealthArray::geom_name(double psi) const
{
  std::stringstream ss;
  ss << "g" ;
  if (psi < 0.1)
  { ss << "0";
    if (psi < 0.01)
    { ss << "0";
      if (psi < 0.001)
      { ss << "0";
        if (psi < 0.0001)
	{ ss << "0";
	  if (psi < 0.00001)
	    ss << "0";
	}
      }
    }
  }
  ss << floor(100000*psi);
  return ss.str();
}



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
  { for (int i=0; i<(int)mRejectPositions.size(); ++i)
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
  else if (0.02 < wealth) return 0.005;
  else return 0.001;
}
