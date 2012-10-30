#include "wealth.h"

#include <iostream>

#include <ctime>
#include <Eigen/Core>
 

// output pair ... why is this wrong?
/*
  inline
  std::ostream&
  std::operator<<(std::ostream& os, std::pair<double,double> const& p)
  { os << "<" << p.first << "," << p.second << ">";
  return os;
  }
*/


int  main()
{
  const int univStart (1);

      
  if (false)
  { double   W0  = 0.05;
    double omega = 0.05;
    double scale =  0.5;   // wealth array gets very large as the scale increases (about 3000 for scale=4, w=0.5; 22000 for scale=0.5,w=0.05)
    std::cout << "\n\nTEST: Scaled wealth (" << scale << ") starting from wealth W0=" << W0
	      << " with universal bids throughout and omega " << omega << std::endl;
    ScaledUniversalDist u(scale);
    int nRounds = 100;
    WealthArray wealth(W0, omega, nRounds, u);
    int iZ (wealth.zero_index());
    std::cout << "   Wealth at position iZero=" << iZ << " is " << wealth.wealth(iZ) << " with next bid to be " << wealth.bid(iZ) << std::endl;
    // std::cout << "    : " << wealth << std::endl;
    for(int i=0; i<10; ++i)
      std::cout << "W[" << i << "] = " << wealth[i] << "  with bid " << wealth.bid(i) << std::endl;
  }

  if (false)
  {
    std::cout << "\n\nTEST: Test bidding from wealth array." << std::endl;
    double omega (  0.05);
    int  nRounds ( 20   );
    int    iZero ( 10   );
    
    Distribution *p;
    UniversalDist univ(univStart);
    p = &univ;
    WealthArray uWealth(omega, iZero, nRounds, *p);
    GeometricDist geo(0.005);
    p = &geo;
    WealthArray gWealth(omega, iZero, nRounds, *p);
    std::cout << "TEST: wealth array  \n" << uWealth << std::endl;
    std::cout << "TEST: wealth array  \n" << gWealth << std::endl;
    std::cout << "TEST: high wealth bids, then those starting from iZero" << std::endl;
    for(int k=0; k<5; ++k)
      std::cout << "  bid at k " << k << " geo bid=" <<  gWealth.bid(k) << " out of " << gWealth[k]
		<< "    universal bid=" << uWealth.bid(k) << " out of " << uWealth[k] << std::endl;
    for (int r=0; r < nRounds; ++r) std::cout << "round=" << r+1
					      << " geo bid=" <<  gWealth.bid(iZero+r) << " out of " << gWealth[iZero+r]
					      << "    universal bid=" << uWealth.bid(iZero+r) << " out of " << uWealth[iZero+r] << std::endl;
  }

  
  if (false)
  {
    std::cout << "\n\nTEST: Test bracketing search in wealth." << std::endl;
    double omega (  0.05);
    int  nRounds ( 50   );
    int    iZero ( 10   );
    UniversalDist univ(1);

    std::vector<int> ii = { 3, 6, 10, 15, 25};
    WealthArray uWealth(omega, iZero, nRounds, univ);
    for (unsigned int j=0; j<ii.size()-1; ++j)
    { int i = ii[j];
      double bid = uWealth.bid(i);
      std::pair<int,double>  kk (uWealth.wealth_position(i));
      double value =  uWealth[kk.first] * kk.second + uWealth[kk.first+1] * (1-kk.second);
      std::cout << "TEST:  increment W[" << i << "]= " << uWealth[i] << " by " << 0.05-bid << " to " << 0.05+uWealth[i]-bid << " bracketed by "
		<< uWealth[kk.first]   << " * (" << kk.second << ")  +  "
		<< uWealth[kk.first+1] << " * (" << 1-kk.second << ") = " << value << std::endl;
    }
  }
  
  if (false)
  { std::cout << "\n\n Test extremes in geometric wealth table for underflows" << std::endl;
    
    double omega ( 0.05 );
    int    iZero (   5  );
    int    steps (  50  );
    double psi(0.01);
      
    UniversalDist univ(univStart);
    GeometricDist geo(psi);
    UniformDist uni(1+iZero);  // add one to cover 0th position
 
    //    WealthArray gWealth(" Geom ", omega, iZero, geo );   // not numerically stable for long trials
    std::cout << "TEST: init universal wealths \n" ;  WealthArray uWealth(omega, iZero, steps, univ);
    std::cout << "TEST: init geometric wealths \n" ;  WealthArray uniformWealth(omega, iZero, steps, uni);
    std::cout << "TEST: init uniform   wealths \n" ;  WealthArray gWealth(omega, iZero, steps, psi );  // better geometric
    
    std::cout << "TEST: geometric name for psi=" << psi << " is " << gWealth.name() << std::endl;

    // wealth at 0 should be 1 step from zero for uniform
    std::cout << "TEST: wealth at  0 is (univ) " << uWealth[ 0]      << "   (geo) " << gWealth[ 0]      << "   (unif) " << uniformWealth[0]     << std::endl;
    std::cout << "TEST: wealth at  1 is        " << uWealth[ 1]      << "         " << gWealth[ 1]      << "          " << uniformWealth[1]     << std::endl;
    std::cout << "TEST: wealth at iZero is     " << uWealth[ iZero ] << "         " << gWealth[ iZero ] << "          " << uniformWealth[iZero] << std::endl;
    std::cout << "TEST: bid(0)                 " << uWealth.bid(0)   << "         " << gWealth.bid(0)   << "          " << uniformWealth.bid(0) << std::endl;

    std::cout << "TEST: Bid comparisons, geometric(0.01) and universal...\n";
    for (int j=1; j<10; ++j)
      std::cout << "[" << j << "]  " << gWealth.bid(iZero-j) << "   " << uWealth.bid(iZero-j) << "   " << uniformWealth.bid(iZero-j) << std::endl;
    
    std::cout << "TEST: universal wealth array  \n" << uWealth << std::endl;
    std::cout << "TEST: geometric wealth array  \n" << gWealth << std::endl;
    std::cout << "TEST:   uniform wealth array  \n" << uniformWealth << std::endl;
  } 

  if (true)
  { std::cout << "\n\nTEST: Testing dual wealth array." << std::endl;
    DualWealthArray wealth("test", 0.5, 0.5, UniversalBidder(2.0), 20);
    std::cout << "TEST: Zero index is at " << wealth.zero_index() << std::endl;
    wealth.write_to(std::cout);
    std::cout << "\n\nTEST: Finished test of dual wealth array" << std::endl;
  }
  
  return 0;
}
