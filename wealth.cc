#include "wealth.h"
#include "line_search.Template.h"

#include <utility>
#include <math.h>
#include <sstream>

static const std::string messageTag ("WLTH: ");

/////////////////////////////////////////  Distributions  ////////////////////////////////////////

double
GeometricDist::operator() (int k) const
{
  double p=1.0;
  for (int i=0; i<k; ++i) p *= mOneMinusPsi;
  return p * mPsi;
}


std::string
GeometricDist::identifier() const
{
  std::stringstream ss;
  ss << "g" << mPsi;
  return ss.str();
}


double
UniformDist::operator() (int k) const
{
  assert ((0 <= k) && (k <= mLimit));
  return mP;
}


std::string
UniformDist::identifier() const
{
  std::stringstream ss;
  ss << "uniform(" << mLimit << ")";
  return ss.str();
}



std::string
UniversalDist::identifier() const
{
  std::stringstream ss;
  ss << "univ(" << mStart << ")";
  return ss.str();
}

// constants to make universal into density, from MMa
double normalizingConstants[21]={0,3.3877355 , 1.3063666 , 0.8920988 , 0.7186514 , 0.6221371,
				   0.5598396 , 0.51582439, 0.48278679, 0.45689505, 0.4359382,
				   0.4185466 , 0.40382391, 0.39115728, 0.38011245, 0.37037246,
 				   0.36170009, 0.35391396, 0.34687281, 0.34046481, 0.33460018};

double
UniversalDist::operator() (int k) const
{
  double ll = log(k+1+mStart);
  return 1.0/( (k+mStart) * ll * ll * normalizingConstants[mStart]);
}



//     scaled universal     scaled universal     scaled universal     scaled universal     scaled universal     scaled universal

const double ScaledUniversalDist::mSumOfRecipLog = 3.387735532;

std::string
ScaledUniversalDist::identifier() const
{
  std::stringstream ss;
  ss << "scaled_univ(" << mScale <<")";
  return ss.str();
}

double
ScaledUniversalDist::operator() (int k) const
{
  assert (0 < k);
  double ll = log(k+1);
  return mScale/(k * ll * ll);
}


int
ScaledUniversalDist::w0_index(double w0) const
{
  double w (max_wealth());
  int j (1);
  while (w0 < w)
  { w -= this->operator()(j);
    ++j;
  }
  return j;
}


//     uniform to end     uniform to end     uniform to end     uniform to end         

double uniform_to_end (int k, int left)         // equal spread over possible locations
{
  return 1.0/(double)(k + left);
}




//     WealthArray     WealthArray     WealthArray     WealthArray     WealthArray     WealthArray     WealthArray     WealthArray

std::pair<int, double>
WealthArray::find_wealth_position (int k0, double increase)  const // k0 is current position denoting current wealth
{
  double target = mWealth[k0] + increase;      // 'wealth' is 'new wealth' > 'current wealth'
  int k1 (mSize-1);                            // W[k0] <= W[k1]
  if (mWealth[k1] < target)                    // outside wealth range, as if bid > payoff
    target = mWealth[k1];
  while (k0+1 < k1)                            // bracket between k0 and k1
  { int kk = floor((k0+k1)/2);
    // std::cout << "DEBUG: find_wealth_position   W[" << k0 << "]="
    //           << mWealth[k0] << " W[" << kk << "]=" << mWealth[kk] << "  W[" << k1 << "]=" << mWealth[k1] << std::endl;
    if (target < mWealth[kk])
    { k1 = kk; }
    else
    { k0 = kk; }
  }
  if (k0<k1)  // inside range
  { double p ( (target - mWealth[k0]) / (mWealth[k1] - mWealth[k0]) );
    if (p < 0) std::cerr << messageTag << "*** Error ***  Wealth position is " << k0 << " " << p << std::endl;
    return std::make_pair(k0,1-p);
  }
  else
    return std::make_pair(k0,1);
}


void
WealthArray::initialize_array_using_pdf(Distribution const& p)
{
  // std::cout << "WARRAY: Building dyn array with " << mSize << " steps and wealth " << mOmega << " @ " << mZeroIndex << std::endl;
  assert((0 < mZeroIndex) && (mZeroIndex < mSize-1));
  DynamicArray<double> da(0,mSize-1);
  da.assign(mZeroIndex,mOmega);
  for(int i=mZeroIndex-1; 0 <= i; --i)
  { // std::cout << " da[i="<<i<<"] = (da[i+1="<<i+1<<"]="<< da[i+1]<<")-("<<mOmega<<")*(p["<< mZeroIndex-i<<"]="<< p(mZeroIndex-i)<<")\n";
    double bid (mOmega * p(mZeroIndex-i-1));    // Note... error would be: mZeroIndex-i 'banks' some wealth
    da.assign(i, da[i+1] - bid);
  }
  mWealth=da;
  fill_array_top();
  init_positions();
}


void
WealthArray::initialize_array_using_func(ScaledUniversalDist const& f)
{
  assert((0 < mZeroIndex) && (mZeroIndex < mSize-1));
  DynamicArray<double> da(0,mSize-1);
  da.assign(mSize-1, f.max_wealth());
  for(int i=1, j=mSize-2; i<mSize; ++i,--j)
    da.assign(j, da[j+1]-f(i));
  mWealth = da;
  init_positions();
}


void
WealthArray::initialize_geometric_array(double psi)
{
  assert((0 < mZeroIndex) && (mZeroIndex < mSize-1));
  DynamicArray<double> da(0,mSize-1);
  da.assign(mZeroIndex,mOmega);
  for(int i=mZeroIndex-1; 0 <= i; --i)
  { double bid (da[i+1]*psi);
    da.assign(i, da[i+1] - bid ); 
  }
  mWealth=da;
  fill_array_top();
  init_positions();
}


void
WealthArray::fill_array_top()
{ if (mPadding > 2)                   // Add padding to accumulate wealth above omega by incrementing omega over padding steps
  { double w (0.5);                   // allow to grow this much
    int    k (mPadding-2) ;           // over this many steps
    // geometric sum
    double b (mWealth[mZeroIndex]-mWealth[mZeroIndex-1]);   // incrementing initial bid
    double m (Line_Search::Bisection(0.00001,std::make_pair(1.000001,3))
	      ([&w,&k,&b](double x){ double xk(x); for(int j=1;j<k;++j) xk *= x; return x*(1.0-xk)/(1-x) - w/b;}));
    if (m < 1)
    { m = 1.0;
      std::cerr << messageTag << " *** Error ***  Wealth array cannot initialize upper wealth for inputs. Setting m = 1." << std::endl;
      std::cout << "            w=" << w << "    k=" << k << "   b=" << b << std::endl;
    }
    for(int i=mZeroIndex+1; i < mSize-1; ++i)
    { b *= m;
      mWealth.assign(i, mWealth[i-1] + b);
    }
  }
  // last increment must be omega
  mWealth.assign(mSize-1, mWealth[mSize-2] + mOmega);
}


void
WealthArray::init_positions ()
{
  // lock in indexing for finding new positions since the increment is known in advance
  mPositions.push_back( std::make_pair(0,0) ) ;
  for(int j = 1; j<mSize-1; ++j)
  { double increase (mOmega - bid(j));
    if (increase < 0)
    { std::cerr << messageTag << "*Warning*  Wealth implies certain loss because bid " << bid(j) << " exceeds payoff " << mOmega << ".  Will stay at current position rather than loss." << std::endl;
      increase = 0;
    }
    mPositions.push_back( find_wealth_position(j,increase) );
  }
}


void
WealthArray::print_to (std::ostream& os) const
{
  os << "Wealth array " << mName << "  has wealth "
     << mWealth[mZeroIndex] << " at iZero=" << mZeroIndex
     << " with wealth vector : \n" << mWealth;
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

