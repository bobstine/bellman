#include "bellman.h"
#include "random.h"


#include <iostream>
#include <algorithm> // generate


inline std::ostream&
operator<<(std::ostream &output, const std::pair<double,double> pair)
{ output << " < " << pair.first << " , " << pair.second << " > "; return output; }



// Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     
int  main()
{
  int    nRounds (100   );
  double omega   (  0.5 );
  double alpha   (  0.05);
  double scale   (  2.0 );
  double angle   (210   );

  // set up utility and wealth
  WealthArray bidderWealth(omega, omega, nRounds, ScaledUniversalDist(scale));
  RiskVectorUtility utility(angle, alpha);

  // generate means
  RandomGenerator randu(23984);
  const double pNull  (0.6);
  const double signal (2.5);
  std::vector<double> means (nRounds);

  const int nReps (50);
  for (int rep = 0; rep < nReps; ++rep)
  { std::generate(means.begin(), means.end(), [&randu,pNull,signal]()->double { if(randu.uniform() < pNull) return 0.0; else return signal; });
    std::pair<double,double> result (find_process_risk (utility, bidderWealth, means));
    std::cout << result.first << " " << result.second << std::endl;
  }
  return 0;
}

