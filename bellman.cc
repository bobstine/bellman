#include "bellman.h"
#include "random.h"
#include "line_search.Template.h"
#include "wealth.h"
#include "eigen_utils.h"


#include <math.h>
#include <functional>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ios>

// #define SHOW_PROGRESS

const std::string messageTag = "BELL: ";

int
imin(int a, int b)
{ if (a < b) return a; else return b; }


Line_Search::GoldenSection
make_search_engine(void)
{
  const int                      maxIterations   (200);   
  const double                   tolerance       (0.0001);
  const double                   initialGridSize (0.5);
  const std::pair<double,double> searchInterval  (std::make_pair(0.05,10.0));
  return Line_Search::GoldenSection(tolerance, searchInterval, initialGridSize, maxIterations);
}


//  Find the risk associated with a Bayesian spike model; this is the code that generates the paths within feasible set

std::pair<double,double>
find_process_risk (int nRounds, double pZero, double mu, VectorUtility & utility, DualWealthArray const& bidderWealth)
{
  const int nColumns (bidderWealth.number_wealth_positions());   
  // pad arrays since need room to collect bid; initialize to zero; extra row for starting bottom up recursion
  Matrix oracleMat = Matrix::Zero(nRounds+1, nColumns+1);
  Matrix bidderMat = Matrix::Zero(nRounds+1, nColumns+1);
  // store intermediates in rectangular array; fill from bottom up
  for (int row = nRounds-1; row > -1; --row)
  { for (int k=0; k<nColumns; ++k)
    { double bid (bidderWealth.bid(k));
      std::pair<int,double> rejectPos (bidderWealth.reject_position(k));          // where to go if reject (col, prob)
      std::pair<int,double>    bidPos (bidderWealth.bid_position(k));             // where to go if dont reject
      double bidderIfReject  = bidderMat(row+1,rejectPos.first)*rejectPos.second + bidderMat(row+1,rejectPos.first+1)*(1-rejectPos.second);
      double oracleIfReject  = oracleMat(row+1,rejectPos.first)*rejectPos.second + oracleMat(row+1,rejectPos.first+1)*(1-rejectPos.second);
      double bidderIfDont    = bidderMat(row+1,   bidPos.first)*   bidPos.second + bidderMat(row+1,   bidPos.first+1)*(1-   bidPos.second);
      double oracleIfDont    = oracleMat(row+1,   bidPos.first)*   bidPos.second + oracleMat(row+1,   bidPos.first+1)*(1-   bidPos.second);
      utility.set_constants(bid, 0.0, 0.0);
      bidderMat (row,k) = pZero * utility.bidder_utility(0, bidderIfReject, bidderIfDont)
	                  + (1-pZero) * utility.bidder_utility(mu, bidderIfReject, bidderIfDont);
      oracleMat (row,k) = pZero * utility.oracle_utility(0, oracleIfReject, oracleIfDont)
 	                  + (1-pZero) * utility.oracle_utility(mu, oracleIfReject, oracleIfDont);
    }
  }
  int iZero (bidderWealth.zero_index());
  return std::make_pair(oracleMat(0,iZero), bidderMat(0,iZero));
}
  

//    solve_bellman_utility dual_wealth     solve_bellman_utility dual_wealth     solve_bellman_utility dual_wealth     solve_bellman_utility dual_wealth
//
//        * Only one is considered to be constrained.  Oracle is automatically unconstrained.
//        * Version uses wealth sliced on y-axis, Lebesgue style
//

void
solve_bellman_utility  (int nRounds, VectorUtility &utility, DualWealthArray const& wealth, bool writeDetails)
{
  const int nColumns (wealth.number_wealth_positions());   
  auto search = make_search_engine();
  // pad arrays since need room to collect bid; initialize to zero; extra row for starting bottom up recursion
  Matrix utilityMat= Matrix::Zero(nRounds+1, nColumns+1);  // +1 for initializing
  Matrix oracleMat = Matrix::Zero(nRounds+1, nColumns+1);
  Matrix bidderMat = Matrix::Zero(nRounds+1, nColumns+1);
  // save for recreating mean stochastic process
  Matrix meanMat       = Matrix::Zero(nRounds  , nColumns);
  Matrix rejectProbMat = Matrix::Zero(nRounds  , nColumns);    // rejection prob
  // store intermediates in rectangular array; fill from bottom up (don't get trapezoid with dual wealth)
  for (int row = nRounds-1; row > -1; --row)
  { for (int k=0; k<nColumns; ++k)
    { double bid (wealth.bid(k));
      std::pair<int,double> rejectPos (wealth.reject_position(k));                 // where to go if reject (col, prob)
      std::pair<int,double> bidPos    (wealth.bid_position(k));                    // did not reject
      double utilityIfReject = utilityMat(row+1,rejectPos.first)*rejectPos.second + utilityMat(row+1,rejectPos.first+1)*(1-rejectPos.second);
      double utilityIfBid    = utilityMat(row+1,   bidPos.first)*   bidPos.second + utilityMat(row+1,   bidPos.first+1)*(1-   bidPos.second);
      double bidderIfReject  =  bidderMat(row+1,rejectPos.first)*rejectPos.second +  bidderMat(row+1,rejectPos.first+1)*(1-rejectPos.second);
      double bidderIfBid     =  bidderMat(row+1,   bidPos.first)*   bidPos.second +  bidderMat(row+1,   bidPos.first+1)*(1-   bidPos.second);
      double oracleIfReject  =  oracleMat(row+1,rejectPos.first)*rejectPos.second +  oracleMat(row+1,rejectPos.first+1)*(1-rejectPos.second);
      double oracleIfBid     =  oracleMat(row+1,   bidPos.first)*   bidPos.second +  oracleMat(row+1,   bidPos.first+1)*(1-   bidPos.second);
      utility.set_constants(bid, utilityIfReject, utilityIfBid);          
      std::pair<double,double> maxPair (search.find_maximum(utility));            // mean and maximal utility
      double utilAtMuEqualZero (utility(0.0));
      if (maxPair.second < utilAtMuEqualZero)
	maxPair = std::make_pair(0.0,utilAtMuEqualZero);
      meanMat       (row,k) = maxPair.first;
      rejectProbMat (row,k) = reject_prob(meanMat(row,k), (bid < 0.99) ? bid : 0.99); // insure prob less than 1
      utilityMat    (row,k) = maxPair.second;
      bidderMat     (row,k) = utility.bidder_utility(maxPair.first, bidderIfReject, bidderIfBid);
      oracleMat     (row,k) = utility.oracle_utility(maxPair.first, oracleIfReject, oracleIfBid);
    }
  }
  // write solution (without boundary row) to file
  if(writeDetails)
  { std::ostringstream ss;
    int angle (trunc(utility.angle()));
    ss << "sim_details/dual_bellman.a" << angle << ".n" << nRounds << ".w" << round(100*wealth.initial_wealth())
       << ".o" << round(100*wealth.omega()) << ".al" << round(100*utility.alpha()) << ".";
    write_matrix_to_file(ss.str() + "utility",    utilityMat.topLeftCorner(nRounds+1,    utilityMat.cols()-1));  // omit boundary row, col
    write_matrix_to_file(ss.str() + "oracle" ,     oracleMat.topLeftCorner(nRounds+1,     oracleMat.cols()-1));
    write_matrix_to_file(ss.str() + "bidder" ,     bidderMat.topLeftCorner(nRounds+1,     bidderMat.cols()-1));
    write_matrix_to_file(ss.str() + "mean"   ,       meanMat.topLeftCorner(nRounds  ,       meanMat.cols()));
    write_matrix_to_file(ss.str() + "prob"   , rejectProbMat.topLeftCorner(nRounds  , rejectProbMat.cols()));
    { std::ios_base::openmode mode = std::ios_base::trunc;
      std::ofstream output (ss.str() + "wealth", mode);
      wealth.write_to(output, true);  // true = as lines
      output.close();
    }
  }
  // locate starting position in array
  int iZero = wealth.zero_index();
  std::cout << std::setprecision(8)
	    << utility.angle() << " " << wealth.omega() << "   " << nRounds   << "   " 
	    << utilityMat(0,iZero) << " " << oracleMat(0,iZero) << " " << bidderMat(0,iZero) << std::endl;
}


//    solve_bellman_utility matrix      solve_bellman_utility matrix     solve_bellman_utility matrix     solve_bellman_utility matrix
//
//  Dual randomized wealth version
//

void
solve_bellman_utility  (int nRounds, MatrixUtility &utility, DualWealthArray const& rowWealth,  DualWealthArray const& colWealth, std::string config, bool writeTable)
{
  const int nRows (1+rowWealth.number_wealth_positions());                // extra 1 for padding; allow 0 * 0
  const int nCols (1+colWealth.number_wealth_positions());
  const std::pair<int,int> zeroIndex(std::make_pair(rowWealth.zero_index(), colWealth.zero_index()));
  std::clog << messageTag <<  "Zero indices are " << zeroIndex.first << " & " << zeroIndex.second << std::endl;
  // code flips between these on read and write using use0
  bool use0 = true;
  Matrix  utilityMat0= Matrix::Zero (nRows, nCols);                       // initialize with zero for risk, 1 for rejection
  Matrix  utilityMat1= Matrix::Zero (nRows, nCols);
  Matrix* pUtilitySrc(&utilityMat0), *pUtilityDest(&utilityMat1);
  Matrix  rowMat0    = Matrix::Zero (nRows, nCols);
  Matrix  rowMat1    = Matrix::Zero (nRows, nCols);
  Matrix* pRowSrc(&rowMat0), *pRowDest(&rowMat1);
  Matrix  colMat0    = Matrix::Zero (nRows, nCols);
  Matrix  colMat1    = Matrix::Zero (nRows, nCols);
  Matrix* pColSrc(&colMat0 ), *pColDest(&colMat1);
  // iteration vars
  std::pair<double,double> maxPair;                                       // x,f(x)
  std::pair<double,double> bestMeanInterval = std::make_pair(10,0);
  auto search = make_search_engine();
  for (int round = nRounds; 0 < round; --round)
  { if (use0)   // flip progress arrays
    { pUtilitySrc = &utilityMat0;    pRowSrc  = &rowMat0;   pColSrc  = &colMat0;
      pUtilityDest= &utilityMat1;    pRowDest = &rowMat1;   pColDest = &colMat1;
    } else
    { pUtilitySrc = &utilityMat1;    pRowSrc  = &rowMat1;   pColSrc  = &colMat1;
      pUtilityDest= &utilityMat0;    pRowDest = &rowMat0;   pColDest = &colMat0;
    }
    use0 = !use0;
    for (int r=0; r<nRows-1; ++r)                                         //  padding... allows zero weight on zero value without if/else
    { double rowBid = rowWealth.bid(r);
      std::pair<int, double> rowBidPos    = rowWealth.bid_position(r);    // if does not reject
      std::pair<int, double> rowRejectPos = rowWealth.reject_position(r); // if rejects
      for (int c=0; c<nCols-1; ++c) 
      { double colBid = colWealth.bid(c);
	std::pair<int, double>  colBidPos   (colWealth.bid_position(c));
	std::pair<int, double> colRejectPos (colWealth.reject_position(c));
	utility.set_constants(rowBid, colBid,
			      reject_value ( rowBidPos  ,  colBidPos  , *pUtilitySrc),   // v00  neither rejects; wealth function is decreasing
			      reject_value ( rowBidPos  , colRejectPos, *pUtilitySrc),   // v01  only column player rejects
			      reject_value (rowRejectPos,  colBidPos  , *pUtilitySrc),   // v10  only row rejects
			      reject_value (rowRejectPos, colRejectPos, *pUtilitySrc));  // v11  both reject
	maxPair = search.find_maximum(utility);       // returns opt, f(opt)
	// monitor range of optimal means
	if(maxPair.first < bestMeanInterval.first)
	  bestMeanInterval.first = maxPair.first;
	else if (maxPair.first > bestMeanInterval.second)
	  bestMeanInterval.second = maxPair.first;	
	double utilAtMuEqualZero = utility(0.0);
	if (maxPair.second < utilAtMuEqualZero)
	  maxPair = std::make_pair(0.0,utilAtMuEqualZero);
	(*pUtilityDest)(r,c) = maxPair.second;
	(* pRowDest)(r,c) = utility.row_utility(maxPair.first,                                            // opt mu
						   reject_value (rowBidPos   ,  colBidPos  , *pRowSrc),   // v00  neither rejects
						   reject_value (rowBidPos   , colRejectPos, *pRowSrc),   // v01  only column player rejects
						   reject_value (rowRejectPos,  colBidPos  , *pRowSrc),   // v10  only row rejects
						   reject_value (rowRejectPos, colRejectPos, *pRowSrc));  // v11  both reject
        (* pColDest)(r,c) = utility.col_utility(maxPair.first,
						   reject_value (rowBidPos   ,  colBidPos  , *pColSrc),   // v00  neither rejects
						   reject_value (rowBidPos   , colRejectPos, *pColSrc),   // v01  only column player rejects
						   reject_value (rowRejectPos,  colBidPos  , *pColSrc),   // v10  only row rejects
						   reject_value (rowRejectPos, colRejectPos, *pColSrc));  // v11  both reject
      }
    }
  }
  std::cout << std::setprecision(6);
  if(writeTable)
  { write_matrix_to_file(config + ".utility", *pUtilityDest);
    write_matrix_to_file(config + ".row" ,  *pRowDest);
    write_matrix_to_file(config + ".col" ,  *pColDest);
  }
  // write summary of configuration and results to stdio
  std::cout << std::setprecision(8)
	    << utility.angle() << " " << 1000*rowWealth.omega()+colWealth.omega() << "   " << nRounds   << "   " 
	    << (*pUtilityDest)(zeroIndex.first, zeroIndex.second) << " "
	    << (*pRowDest    )(zeroIndex.first, zeroIndex.second) << " "
	    << (*pColDest    )(zeroIndex.first, zeroIndex.second) << std::endl;
}


