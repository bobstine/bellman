#include "wealth.h"

#include <iostream>
#include <vector>
#include <ctime>
 

int  main()
{

  const int univStart (1);
  
  if (false)
  {
    std::cout << "\nTest the probability function from wealth.h" << std::endl;
    UniversalDist univ(univStart);
    double total (0.0);
    int count = 100000;
    std::cout << "TEST: initial 20 universal rates (" << univ.identifier() << ")  ";
    for(int k=0; k<20; ++k) std::cout << univ(k) << " "; std::cout << std::endl;
    for(int k=0; k<count; ++k) total += univ(k);
    std::cout << "TEST: Total of universal(*,0.05) for " << count << " terms = " << total << std::endl;

    total = 0.0;
    count = 10000;
    GeometricDist geo(0.005);
    std::cout << "TEST: initial 20 geometric rates (" << geo.identifier() << ")  ";
    for(int k=0; k<20; ++k) std::cout << geo(k) << " "; std::cout << std::endl;
    for(int k=0; k<count; ++k) total += geo(k);
    std::cout << "TEST: Total of geometric for " << count << " terms = " << total << std::endl;

    total = 0.0;
    count = 100;
    UniformDist uni(count);
    std::cout << "TEST: initial 20 uniform rates (" << uni.identifier() << ")  ";
    for(int k=0; k<20; ++k) std::cout << uni(k) << " "; std::cout << std::endl;
    for(int k=0; k<count; ++k) total += uni(k);
    std::cout << "TEST: Total of uniform for " << count << " terms = " << total << std::endl;
  } 

 
  if (true)
  { std::cout << "\n\n--- Test scaled universal wealth function ---" << std::endl;
    std::vector<double> scales = {1,2,4};

    for (unsigned int i=0; i<scales.size(); ++i)
    { double scale = scales[i];
      ScaledUniversalDist u(scale);
      std::cout << "Total wealth for scale " << scale << " is " << u.max_wealth() << std::endl;
      for (int j=1; j<10; ++j)
	std::cout << "    Scale u(" << scale << ")@" << j << " = " << u(j) << std::endl;
      double W0 (0.5);
      std::cout << "    Starting index for initial wealth " << W0 << " is " << u.w0_index(W0) << std::endl;
      W0 = 1.5;
      std::cout << "    Starting index for initial wealth " << W0 << " is " << u.w0_index(W0) << std::endl;
    }
  }

  return 0;
}
