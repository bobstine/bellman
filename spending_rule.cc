#include "spending_rule.h"

#include <math.h>
#include <sstream>

static const std::string messageTag ("SPRL: ");

//     geometric     geometric     geometric     geometric     geometric     geometric

std::string
GeometricRule::identifier() const
{
  std::stringstream ss;
  ss << "Geom(" << mPsi << ")";
  return ss.str();
}

double
GeometricDist::operator() (double w) const
{
  return w * mPsi;
}

//     universal     universal     universal     universal     universal     universal     

std::string
UniversalRule::identifier () const
{
  return "Univ()";
}

static double ln2 =0.69314718055994530942;

			   double
UniversalRule::operator()(double w) const
{
  return w-ln2/log(1.0+exp(ln2/w));
}

