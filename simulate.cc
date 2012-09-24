#include "bellman.h"
#include "utility.h"
#include "random.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include <vector>


RandomGenerator randu(1235);

double
simulate_next_mean()
{
  if (randu.uniform() < 0.5)
    return 2.2;
  else
    return 0.0;
}


void
main (void)
{
  int    nRounds (100);
  double  omega  (0.5);
  double  scale  (1.0);
  double  alpha  (0.05);  // bid of oracle
  
  const int iOmega  (nRounds + 1);
  ScaledUniversalDist univDist = WealthArray(omega, omega, iOmega, ScaledUniversalDist(scale));

  // fill mean vector
  std::vector<double> mu (nRounds);
  for (int i=0; i<nRounds; ++i) mu[i] = simulate_next_mean();

  // find expected oracle risk
  double oracleRisk (0.0);
  for (int t=0; t<nRounds; ++t)
    oracleRisk += risk(mu[t],alpha);

  // simulate bidder risk using nReps
  int nReps (100);
  std::vector<double> bidderRisk(nReps);
  bool bidderRejects (false);
  for (rep = 0; rep < nReps; ++rep)
  { double rBidder = 0;
    double wealth = 0.5;
    int sinceReject = 0;
    for (int t=0; t<nRounds; ++t)
    { bidderRejects = false;
      double beta   (wealth * univDist(sinceReject)) ;
      wealth -= beta;
      double pBeta  (reject_prob(mu,beta ));
      double uni    (randu.uniform());
      if (uni < beta) // reject, estimate by mu
      { rBidder += 1.0;
	wealth += omega;
      }
      else
	rBidder +=mu[i] * mu[i]; 
    }
    bidderRisk[rep] = rBidder;
  }
  std::cout << "Oracle risk " << oracleRisk << "  Bidder risk " << mean(bidderRisk) << std::endl;
}

