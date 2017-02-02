#include "utility.h"

// -------------------------------------------------------------------------------------------------------------

//   RejectMatrixUtility     RejectUtility     RejectUtility     RejectUtility     RejectUtility     RejectUtility     

template<class C>
double
RejectMatrixUtility<C>::operator()(double mu) const
{
  std::pair<double,double>  rprob  (reject_probabilities(mu));
  double rAlpha (rprob.first);
  double rBeta (rprob.second);
  double util (mCriterion(rAlpha,rBeta));
  if (rAlpha > rBeta)
    return  util + mV00 * (1-rAlpha) + mV10 * (rAlpha-rBeta) +  mV11 * rBeta;
  else
    return  util + mV00 * (1- rBeta) + mV01 * (rBeta-rAlpha) +  mV11 * rAlpha;
}

template<class C>
double
RejectMatrixUtility<C>::row_utility (double mu, double v00, double v01, double v10, double v11) const
{
  std::pair<double,double>  rprob  (reject_probabilities(mu));
  double rAlpha (rprob.first);
  double rBeta (rprob.second);
  if (rAlpha > rBeta)
    return  rAlpha + v00 * (1-rAlpha) + v10 * (rAlpha-rBeta) +  v11 * rBeta;
  else
    return  rAlpha + v00 * (1- rBeta) + v01 * (rBeta-rAlpha) +  v11 * rAlpha;
}

template<class C>
double
RejectMatrixUtility<C>::col_utility (double mu, double v00, double v01, double v10, double v11) const
{
  // same recursive structure as oracle_utility, just different recursive values appear in call
  std::pair<double,double>  rprob  (reject_probabilities(mu));
  double rAlpha (rprob.first);
  double rBeta (rprob.second);
  if (rAlpha > rBeta)
    return  rBeta + v00 * (1-rAlpha) + v10 * (rAlpha-rBeta) +  v11 * rBeta;
  else
    return  rBeta + v00 * (1- rBeta) + v01 * (rBeta-rAlpha) +  v11 * rAlpha;
}


template <class C>
void
RejectMatrixUtility<C>::write_details_to_stream(double mu, std::ostream& os) const
{
  std::pair<double, double> rprob  (reject_probabilities(mu));
  os << "RejectMatrixUtility " << identifier() << " @ mu=" << mu
     << ":  alphaBeta = (" << mAlpha << "," << mBeta << ") "
     << "  reject prob = (" << rprob.first << "," << rprob.second << ") "
     << " P(reject @ mu=" << mu << ",a)=" << reject_prob(mu,mAlpha)
     << " P(reject @ mu=" << mu << ",b)=" << reject_prob(mu,mBeta)
     << " V00=" << mV00 << "  V10=" << mV10 << " V01=" << mV01 << " V11=" << mV11;
}




//    RiskUtility      RiskUtility      RiskUtility      RiskUtility      RiskUtility      RiskUtility      RiskUtility      

template <class C>
double
RiskMatrixUtility<C>::operator()(double mu) const
{
  std::pair<double,double>  rprob  (MatrixUtility::reject_probabilities(mu));
  double rAlpha (rprob.first );
  double rBeta  (rprob.second);
  double util   (mCriterion(risk(mu,MatrixUtility::mAlpha),risk(mu,MatrixUtility::mBeta)));
  if (rAlpha > rBeta)
    return  util + MatrixUtility::mV00 * (1-rAlpha) + MatrixUtility::mV10 * (rAlpha-rBeta) +  MatrixUtility::mV11 * rBeta;
  else
    return  util + MatrixUtility::mV00 * (1- rBeta) + MatrixUtility::mV01 * (rBeta-rAlpha) +  MatrixUtility::mV11 * rAlpha;
}


template<class C>
double
RiskMatrixUtility<C>::row_utility  (double mu, double v00, double v01, double v10, double v11) const
{
  std::pair<double,double>  rprob  (reject_probabilities(mu));
  double rAlpha (rprob.first );
  double rBeta  (rprob.second);
  double util (risk(mu,mAlpha));
  if (rAlpha > rBeta)
    return  util  + v00 * (1-rAlpha) + v10 * (rAlpha-rBeta) +  v11 * rBeta;
  else
    return  util  + v00 * (1- rBeta) + v01 * (rBeta-rAlpha) +  v11 * rAlpha;
}


template <class C>
double
RiskMatrixUtility<C>::col_utility (double mu, double v00, double v01, double v10, double v11) const
{
  std::pair<double,double>  rprob  (reject_probabilities(mu));
  double rAlpha (rprob.first);
  double rBeta (rprob.second);
  double util  (risk(mu,mBeta));
  if (rAlpha > rBeta)
    return  util  + v00 * (1-rAlpha) + v10 * (rAlpha-rBeta) +  v11 * rBeta;
  else
    return  util  + v00 * (1- rBeta) + v01 * (rBeta-rAlpha) +  v11 * rAlpha;
}

template <class C>
void
RiskMatrixUtility<C>::write_details_to_stream(double mu, std::ostream& os) const
{
  std::pair<double, double> rprob  (MatrixUtility::reject_probabilities(mu));
  os << "RiskMatrixUtility " << identifier() << " @ mu=" << mu << "   " << 10*mu 
     << "  rProb = (" << rprob.first << "," << rprob.second << ")   "
     << " P(reject @ mu=" << mu << ",a=" << mAlpha << ")=" << reject_prob(mu,mAlpha)
     << " P(reject @ mu=" << mu << ",b=" << mBeta  << ")=" << reject_prob(mu,mBeta)
     << "  risk(mu,mAlpha)=" << risk(mu,mAlpha) << " risk(mu,mBeta)=" << risk(mu,mBeta);
    //     << " V00=" << mV00 << "  V10=" << mV10 << " V01=" << mV01 << " V11=" << mV11;
}

