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
  std::pair<double,double>  rprob  (MatrixUtility::reject_probabilities(mu));
  double rAlpha (rprob.first );
  double rBeta  (rprob.second);
  double util (risk(mu,MatrixUtility::mAlpha));
  if (rAlpha > rBeta)
    return  util  + v00 * (1-rAlpha) + v10 * (rAlpha-rBeta) +  v11 * rBeta;
  else
    return  util  + v00 * (1- rBeta) + v01 * (rBeta-rAlpha) +  v11 * rAlpha;
}


template <class C>
double
RiskMatrixUtility<C>::col_utility (double mu, double v00, double v01, double v10, double v11) const
{
  std::pair<double,double>  rprob  (MatrixUtility::reject_probabilities(mu));
  double rAlpha (rprob.first);
  double rBeta (rprob.second);
  double util  (risk(mu,MatrixUtility::mBeta));
  if (rAlpha > rBeta)
    return  util  + v00 * (1-rAlpha) + v10 * (rAlpha-rBeta) +  v11 * rBeta;
  else
    return  util  + v00 * (1- rBeta) + v01 * (rBeta-rAlpha) +  v11 * rAlpha;
}



