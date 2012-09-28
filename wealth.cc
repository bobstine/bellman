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
  return mWealth[k]-mWealth[k+1];
}

std::pair<int, double>
WealthArray::find_wealth_position (int k1, double increase)  const // k0 is current position denoting current wealth
{
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
    if (p < 0) std::cerr << messageTag << "*** Error ***  Wealth position is " << k0 << " " << p << std::endl;
    return std::make_pair(k0,1-p);
  }
  else
    return std::make_pair(k0,1);
}


void
WealthArray::init_check() const
{
  std::cout << messageTag << "Initializing wealth array with " << size()-mZeroIndex
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
  // std::cout << messageTag << "Geometric to top, b=" << b << " and growth factor m= " << m << " applied to " << mWealth[mZeroIndex] <<std::endl;
  if (m < 1)
  { m = 1.0;
    std::cerr << messageTag << " *** Error ***  Wealth array cannot initialize upper wealth for inputs. Setting m = 1." << std::endl;
    std::cout << "            w=" << w << "    k=" << k << "   b=" << b << std::endl;
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
  // cache indexing for new positions since the increment omega is fixed
  for(int j = 0; j<size()-1; ++j)
  { double increase = mOmega - bid(j);
    if (increase < 0)
      std::cerr << messageTag << "*Warning*  Wealth implies certain loss because bid " << bid(j) << " exceeds payoff " << mOmega  << std::endl;
    mPositions.push_back( find_wealth_position(j,increase) );
  }
  for (int i=0; i<size()-1; ++i) 
  { std::cout << messageTag << "DEBUG i=" << i << "  $" <<mWealth[i] << "  " << bid(i)
	      << "  [" << mPositions[i].first << "," << mPositions[i].second << "]" << std::endl;
  }
}


void
WealthArray::print_to (std::ostream& os) const
{
  os << "Wealth array " << mName << "  has wealth " << mWealth[mZeroIndex] << " at iZero=" << mZeroIndex
     << " with wealth vector : \n";
  for (int i=0; i<size(); ++i) os << mWealth[i] << " ";
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

