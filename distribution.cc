#include "distribution.h"

#include <utility>
#include <math.h>
#include <assert.h>
#include <sstream>
#include <iostream>

static const std::string messageTag ("DIST: ");

//     geometric     geometric     geometric     geometric     geometric     geometric

std::string
GeometricDist::identifier() const
{
  std::stringstream ss;
  ss << "g" << mPsi;
  return ss.str();
}

double
GeometricDist::operator() (int k) const
{
  assert (0 <= k);
  double p = 1/(1-mPsi);
  for (int i=0; i<=k; ++i) p *= mOneMinusPsi;
  return p * mPsi;
}

//     uniform    uniform    uniform    uniform    uniform    uniform    uniform    uniform    

std::string
UniformDist::identifier() const
{
  std::stringstream ss;
  ss << "uniform(" << mLimit << ")";
  return ss.str();
}

double
UniformDist::operator() (int k) const
{
  assert ((0 <= k) && (k <= mLimit));
  return mP;
}


//     universal     universal     universal     universal     universal     universal     

// constants to make universal into density depending on starting index (from MMa)
double normalizingConstants[21]={0,3.3877355 , 1.3063666 , 0.8920988 , 0.7186514 , 0.6221371,
				   0.5598396 , 0.51582439, 0.48278679, 0.45689505, 0.4359382,
				   0.4185466 , 0.40382391, 0.39115728, 0.38011245, 0.37037246,
 				   0.36170009, 0.35391396, 0.34687281, 0.34046481, 0.33460018};

std::string
UniversalDist::identifier() const
{
  std::stringstream ss;
  ss << "univ(" << mStart << ")";
  return ss.str();
}

double
UniversalDist::operator() (int k) const
{
  assert (0 <= k);
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
  assert (0 <= k);
  ++k;
  double ll = log(k+1);
  return mScale/(k * ll * ll);
}


int
ScaledUniversalDist::w0_index(double w0) const
{
  const int sizeLimit (100000);
  double w (max_wealth());
  int j (0);
  while ((w0 < w) && (j < sizeLimit))
  { w -= this->operator()(j);
    ++j;
  }
  if(sizeLimit == j)
    std::clog << messageTag << "*** Error *** Scaled distribution requires more than " << sizeLimit << " positions." << std::endl;
  return j;
}


//     uniform to end     uniform to end     uniform to end     uniform to end         

double uniform_to_end (int k, int left)         // equal spread over possible locations
{
  return 1.0/(double)(k + left);
}



