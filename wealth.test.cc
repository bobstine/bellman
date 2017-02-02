#include "wealth.Template.h"

#include <iostream>

#include <ctime>
#include <Eigen/Core>
 


int  main()
{
      
  if (true)
  { std::cout << "\n\nTEST: Testing dual wealth array, W0=omega=0.25." << std::endl;
    double omega = 0.25;
    int rounds   = 1000;
    double maxWealth = 5.0;
    double psi = 0.001;
    DualWealthArray bonferroni(0.001);
    std::cout << "\nTEST: Bonferroni wealth array: \n ";
    bonferroni.write_to(std::cout);
    std::cout << "\nTEST: Universal wealth array: \n";
    DualWealthArray univWealth("univ", maxWealth, omega, omega, UniversalRule(), rounds);
    std::cout << "TEST: Zero index is at " << univWealth.zero_index() << std::endl;
    univWealth.write_to(std::cout);
    std::cout << "TEST: Geometric wealth array: \n";
    DualWealthArray geoWealth("geo", maxWealth, omega, omega, GeometricRule(psi), rounds);
    std::cout << "TEST: Zero index is at " << geoWealth.zero_index() << std::endl;
    geoWealth.write_to(std::cout);
    std::cout << "------------- " << std::endl << "TEST: Finished test of dual wealth array" << std::endl;
  }

  
  if (false)
  { std::cout << "\n\nTEST: Testing dual wealth array with initial value W0=0.7 != omega=0.05 (compare above)." << std::endl;
    double w0 = 0.70;
    double omega = 0.05;
    int rounds   = 20;
    double maxWealth = 5.;
    double rate = 0.1;
    DualWealthArray geoWealth("geo", maxWealth, w0, omega, GeometricRule(rate), rounds);
    std::cout << "TEST: Zero index is at " << geoWealth.zero_index() << std::endl;
    geoWealth.write_to(std::cout);
    std::cout << "------------- " << std::endl << "TEST: Finished test of dual wealth array with W0 != omega" << std::endl;
  }

  return 0;
}
