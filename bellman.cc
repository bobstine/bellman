#include "bellman.h"
#include "random.h"
#include "line_search.Template.h"
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
  const double                   initialGrid     (0.5);
  const std::pair<double,double> searchInterval  (std::make_pair(0.05,10.0));
  return Line_Search::GoldenSection(tolerance, searchInterval, initialGrid, maxIterations);
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
  

std::pair<double,double>
find_process_risk (int nRounds, double pZero, double mu, VectorUtility & utility, WealthArray const& bidderWealth)
{
  const int nColumns (bidderWealth.number_of_bids());   
  // pad arrays since need room to collect bid; initialize to zero; extra row for starting bottom up recursion
  Matrix oracleMat = Matrix::Zero(nRounds+1, nColumns+1);
  Matrix bidderMat = Matrix::Zero(nRounds+1, nColumns+1);
  // store intermediates in trapezoidal array with shape |\ ; fill from bottom up
  int done = 0;
  for (int row = nRounds-1; row > -1; ++done, --row)
  { for (int k=0; k<nColumns-done; ++k)
    { double bid (bidderWealth.bid(k));
      std::pair<int,double> kp (bidderWealth.wealth_position(k));                 // where to go if reject (col, prob)
      double bidderIfReject  =  bidderMat(row+1,kp.first)*kp.second +  bidderMat(row+1,kp.first+1)*(1-kp.second);
      double oracleIfReject  =  oracleMat(row+1,kp.first)*kp.second +  oracleMat(row+1,kp.first+1)*(1-kp.second);
      utility.set_constants(bid, 0.0, 0.0);
      bidderMat (row,k) = pZero * utility.bidder_utility(0, bidderIfReject, bidderMat(row+1,k+1))
	                  + (1-pZero) * utility.bidder_utility(mu, bidderIfReject, bidderMat(row+1,k+1));
      oracleMat (row,k) = pZero * utility.oracle_utility(0, oracleIfReject, oracleMat(row+1,k+1))
	+ (1-pZero) * utility.oracle_utility(mu, oracleIfReject, oracleMat(row+1,k+1));
    }
  }
  int iZero = bidderWealth.number_of_bids() - nRounds;
  return std::make_pair(oracleMat(0,iZero), bidderMat(0,iZero));
}
  

//
//    Unconstrained   Unconstrained   Unconstrained   Unconstrained   Unconstrained   Unconstrained
//
//    solve_bellman_utility dual_wealth     solve_bellman_utility dual_wealth     solve_bellman_utility dual_wealth     solve_bellman_utility dual_wealth
//
//       Version uses wealth sliced on y-axis, Lebesgue style
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
    ss << "sim_details/dual_bellman.a" << angle << ".n" << nRounds
       << ".s" << round(10*wealth.scale()) << ".o" << round(100*wealth.omega()) << ".al" << round(100*utility.alpha()) << ".";
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
  std::cout << utility.angle() << " " << wealth.omega() << "   " << nRounds   << "   " 
	    << utilityMat(0,iZero) << " " << oracleMat(0,iZero) << " " << bidderMat(0,iZero) << std::endl;
}



//    solve_bellman_utility     solve_bellman_utility     solve_bellman_utility     solve_bellman_utility     
//
//    This version has the 'Riemann' style bidding function
//

void
solve_bellman_utility  (int nRounds, VectorUtility &utility, WealthArray const& bidderWealth, bool writeDetails)
{
  const int nColumns (bidderWealth.number_of_bids());   
  auto search = make_search_engine();
  // pad arrays since need room to collect bid; initialize to zero; extra row for starting bottom up recursion
  Matrix utilityMat= Matrix::Zero(nRounds+1, nColumns+1);  // +1 for initializing
  Matrix oracleMat = Matrix::Zero(nRounds+1, nColumns+1);
  Matrix bidderMat = Matrix::Zero(nRounds+1, nColumns+1);
  // save for recreating mean stochastic process
  Matrix probMat   = Matrix::Zero(nRounds  , nColumns);    // rejection prob
  Matrix indxMat   = Matrix::Zero(nRounds  , nColumns);    // if reject, where to
  Matrix meanMat   = Matrix::Zero(nRounds  , nColumns);
  // store intermediates in trapezoidal array with shape |\ ; fill from bottom up
  int done = 0;
  for (int row = nRounds-1; row > -1; ++done, --row)
  { for (int k=0; k<nColumns-done; ++k)
    { double bid (bidderWealth.bid(k));
      std::pair<int,double> kp (bidderWealth.wealth_position(k));                 // where to go if reject (col, prob)
      double utilityIfReject = utilityMat(row+1,kp.first)*kp.second + utilityMat(row+1,kp.first+1)*(1-kp.second);
      double bidderIfReject  =  bidderMat(row+1,kp.first)*kp.second +  bidderMat(row+1,kp.first+1)*(1-kp.second);
      double oracleIfReject  =  oracleMat(row+1,kp.first)*kp.second +  oracleMat(row+1,kp.first+1)*(1-kp.second);
      utility.set_constants(bid, utilityIfReject, utilityMat(row+1,k+1));         // last is util if not reject
      std::pair<double,double> maxPair (search.find_maximum(utility));            // mean and maximal utility
      double utilAtMuEqualZero (utility(0.0));
      if (maxPair.second < utilAtMuEqualZero)
	maxPair = std::make_pair(0.0,utilAtMuEqualZero);
      meanMat   (row,k) = maxPair.first;
      probMat   (row,k) = reject_prob(meanMat(row,k), (bid < 0.99) ? bid : 0.99); // insure prob less than 1
      indxMat   (row,k) = kp.first;                                               // ignore averaging in rejection destination
      utilityMat(row,k) = maxPair.second;
      bidderMat (row,k) = utility.bidder_utility(maxPair.first, bidderIfReject, bidderMat(row+1,k+1));
      oracleMat (row,k) = utility.oracle_utility(maxPair.first, oracleIfReject, oracleMat(row+1,k+1));
      if (false) std::cout << "     @ " << k << " " << row << " k_r= " << kp.first << " { (" << utilityIfReject << "="
			   << utilityMat(row+1,kp.first  ) << "*" <<   kp.second << " + "
			   << utilityMat(row+1,kp.first+1) << "*" << 1-kp.second  << "), " << utilityMat(row+1,k+1)
			   << "}  with bid " << bidderWealth.bid(k) << "    max  : " << maxPair.second << " @ " << maxPair.first << std::endl;
    }
  }
  // write solution (without boundary row) to file
  if(writeDetails)
  { std::ostringstream ss;
    int angle (trunc(utility.angle()));
    ss << "sim_details/bellman.a" << angle << ".n" << nRounds
       << ".o" << round(100*bidderWealth.omega()) << ".al" << round(100*utility.alpha()) << ".";
    write_matrix_to_file(ss.str() + "utility", utilityMat.topLeftCorner(nRounds+1, utilityMat.cols()-1));  // omit boundary row, col
    write_matrix_to_file(ss.str() + "oracle" ,  oracleMat.topLeftCorner(nRounds+1, oracleMat.cols()-1));
    write_matrix_to_file(ss.str() + "bidder" ,  bidderMat.topLeftCorner(nRounds+1, bidderMat.cols()-1));
    write_matrix_to_file(ss.str() + "mean"   ,    meanMat.topLeftCorner(nRounds  , meanMat.cols()));
    write_matrix_to_file(ss.str() + "prob"   ,    probMat.topLeftCorner(nRounds  , probMat.cols()));
    write_matrix_to_file(ss.str() + "indx"   ,    indxMat.topLeftCorner(nRounds  , indxMat.cols()));
    { std::ios_base::openmode mode = std::ios_base::trunc;
      std::ofstream output (ss.str() + "wealth", mode);
      bidderWealth.print_to(output);
      bidderWealth.write_to(output);
      output.close();
    }
  }
  // locate starting position in |\ array
  int iZero = bidderWealth.number_of_bids() - nRounds;
  std::cout << utility.angle() << " " << bidderWealth.omega() << "   " << nRounds   << "   "
	    << utilityMat(0,iZero) << " " << oracleMat(0,iZero) << " " << bidderMat(0,iZero) << std::endl;
}



//    solve_bellman_utility matrix      solve_bellman_utility matrix     solve_bellman_utility matrix     solve_bellman_utility matrix
//
//  Dual randomized wealth version
//

void
solve_bellman_utility  (int nRounds, MatrixUtility &utility, WealthArray const& rowWealth,  DualWealthArray const& colWealth, bool writeDetails)
{
  const int nRows (rowWealth.number_of_bids());   
  const int nCols (colWealth.number_wealth_positions());   
  // code flips between these on read and write using use0
  std::cout << messageTag <<  " Allocating space " << nRows << "   &   " << nCols << std::endl;
  bool use0 = true;
  Matrix  utilityMat0= Matrix::Zero (nRows, nCols);       // initialize with zero for risk, 1 for rejection
  Matrix  utilityMat1= Matrix::Zero (nRows, nCols);
  Matrix* pUtilitySrc (&utilityMat0), * pUtilityDest (&utilityMat1);
  Matrix  rowMat0    = Matrix::Zero (nRows, nCols);
  Matrix  rowMat1    = Matrix::Zero (nRows, nCols);
  Matrix* pRowSrc (&rowMat0 ), * pRowDest  (&rowMat1);
  Matrix  colMat0    = Matrix::Zero (nRows, nCols);
  Matrix  colMat1    = Matrix::Zero (nRows, nCols);
  Matrix* pColSrc (&colMat0 ), * pColDest  (&colMat1);
  // iteration vars
  std::pair<double,double> maxPair, bestMeanInterval;
  bestMeanInterval = std::make_pair(10,0);
  auto search = make_search_engine();
  for (int round = nRounds; round > 0; --round)
  { std::cout << messageTag << "\n Top of round " << round << std::endl;
    if (use0)   // flip progress arrays
    { pUtilitySrc = &utilityMat0;    pRowSrc  = &rowMat0;   pColSrc  = &colMat0;
      pUtilityDest= &utilityMat1;    pRowDest = &rowMat1;   pColDest = &colMat1;
    } else
    { pUtilitySrc = &utilityMat1;    pRowSrc  = &rowMat1;   pColSrc  = &colMat1;
      pUtilityDest= &utilityMat0;    pRowDest = &rowMat0;   pColDest = &colMat0;
    }
    use0 = !use0;
    for (int r=0; r<nRows-1; ++r) 
    { double rowBid = rowWealth.bid(r);
      std::pair<int, double> rowPos = rowWealth.wealth_position(r); // if rejects
      //      std::cout << messageTag << "Filling row " << r << " and row bid=" << rowBid << " and reject position (prob)" << rowPos.first << " " << rowPos.second << std::endl;
      for (int c=0; c<nCols-1; ++c) 
      { double colBid = colWealth.bid(c);
	std::pair<int, double>  colBidPos   (colWealth.bid_position(c));
	std::pair<int, double> colRejectPos (colWealth.reject_position(c));
	// std::cout << messageTag << "Filling col " << c << " and col bid=" << colBid << " with positions " << colBidPos.first << " " << colRejectPos.first << std::endl;
	utility.set_constants(rowBid, colBid,
			      reject_value ( r+1  ,  colBidPos  , *pUtilitySrc),   // v00  neither rejects; wealth function is decreasing
			      reject_value ( r+1  , colRejectPos, *pUtilitySrc),   // v01  only column player rejects
			      reject_value (rowPos,  colBidPos  , *pUtilitySrc),   // v10  only row rejects
			      reject_value (rowPos, colRejectPos, *pUtilitySrc));  // v11  both reject
	maxPair = search.find_maximum(utility);       // returns opt, f(opt)
	// monitor range of optimal means
	if(maxPair.first < bestMeanInterval.first)
	  bestMeanInterval.first = maxPair.first;
	else if (maxPair.first > bestMeanInterval.second)
	  bestMeanInterval.second = maxPair.first;
	double utilAtMuEqualZero = utility(0.0);
	if (maxPair.second < utilAtMuEqualZero)
	  maxPair = std::make_pair(0.0,utilAtMuEqualZero);
	// std::cout << messageTag << "Utility at (" << r << "," << c << ")= = " << maxPair.first << "  " << maxPair.second << std::endl;
	(*pUtilityDest)(r, c) = maxPair.second;
	(* pRowDest)(r,c) = utility.oracle_utility(maxPair.first,
						   reject_value ( r+1  ,  colBidPos  , *pRowSrc),   // v00  neither rejects
						   reject_value ( r+1  , colRejectPos, *pRowSrc),   // v01  only column player rejects
						   reject_value (rowPos,  colBidPos  , *pRowSrc),   // v10  only row rejects
						   reject_value (rowPos, colRejectPos, *pRowSrc));  // v11  both reject
        (* pColDest)(r,c) = utility.bidder_utility(maxPair.first,
						   reject_value ( r+1  ,  colBidPos  , *pColSrc),   // v00  neither rejects
						   reject_value ( r+1  , colRejectPos, *pColSrc),   // v01  only column player rejects
						   reject_value (rowPos,  colBidPos  , *pColSrc),   // v10  only row rejects
						   reject_value (rowPos, colRejectPos, *pColSrc));  // v11  both reject
	// fill padding column
	(*pUtilityDest)(r,nCols-1) = (* pUtilityDest)(r,nCols-2);
	(*  pRowDest  )(r,nCols-1) = (*    pRowDest )(r,nCols-2);
	(*  pColDest  )(r,nCols-1) = (*    pColDest )(r,nCols-2);
      }
    }
    // fill padding row
    for (int kb=0; kb<nCols; ++kb) 
    { (*pUtilityDest)(nRows-1,nCols-1) = (* pUtilityDest)(nRows-2,nCols-1);
      (*  pRowDest  )(nRows-1,nCols-1) = (*   pRowDest  )(nRows-2,nCols-1);
      (*  pColDest  )(nRows-1,nCols-1) = (*   pColDest  )(nRows-2,nCols-1);
    }
  }
  std::cout << std::setprecision(6);
  if(writeDetails)
  { std::ostringstream ss;
    int angle (utility.angle());
    ss << "runs/bellman2.g" << angle << ".n" << nRounds << ".";
    { write_matrix_to_file(ss.str() + "utility", *pUtilityDest);
      write_matrix_to_file(ss.str() + "row" ,  *pRowDest);
      write_matrix_to_file(ss.str() + "col" ,  *pColDest);
    }
  }
  // write summary of configuration and results to stdio
  std::cout << utility.angle() << " " << rowWealth.omega() << "   " << nRounds   << "   " 
	    << (*pUtilityDest)(nRounds,nRounds) << " " << (*pRowDest)(nRounds,nRounds) << " " << (*pColDest)(nRounds,nRounds) << std::endl;
}


void
solve_bellman_utility  (int nRounds, MatrixUtility & utility, WealthArray const& oracleWealth, WealthArray const& bidderWealth, bool writeDetails)
{
  // initialize: omega location, size includes padding for wealth above omega
  const int iOmega   (nRounds + 1);   
  const int nColumns (bidderWealth.number_of_bids());   
  // pad arrays since need room to collect bid; initialize to zero
  // code flips between these on read and write using use0
  bool use0 = true;
  Matrix utilityMat0= Matrix::Zero (nColumns, nColumns);     // initialize with zero for risk, 1 for rejection
  Matrix utilityMat1= Matrix::Zero (nColumns, nColumns);
  Matrix oracleMat0 = Matrix::Zero (nColumns, nColumns);
  Matrix oracleMat1 = Matrix::Zero (nColumns, nColumns);
  Matrix bidderMat0 = Matrix::Zero (nColumns, nColumns);
  Matrix bidderMat1 = Matrix::Zero (nColumns, nColumns);
  // save two slices of the 'pyramid' of optimal means picked by oracle
  // arrays to save mean chosen by oracle at two fixed wealths identified mIndexA and mIndexB
  // top row A holds bid, second row the wealth for oracle, with specific bid in last col
  // top row B holds information for the bidder
  const int mIndexA = iOmega -  1;
  const int mIndexB = iOmega - imin(50,iOmega/2);  // less wealth
  Matrix meanMatA    = Matrix::Zero (nRounds+2, nColumns);
  Matrix meanMatB    = Matrix::Zero (nRounds+2, nColumns);
  for (int col=0; col<nColumns-1; ++col)
  { meanMatA(0,col) = oracleWealth.bid(col);
    meanMatA(1,col) = oracleWealth[col];
    meanMatB(0,col) = bidderWealth.bid(col);
    meanMatB(1,col) = bidderWealth[col];
  }
  meanMatA(0,nColumns-1) = oracleWealth.bid(mIndexA);
  meanMatB(0,nColumns-1) = oracleWealth.bid(mIndexB);
  // alternate between reading and writing these matrices
  Matrix* pUtilitySrc (&utilityMat0), * pUtilityDest (&utilityMat1);
  Matrix* pOracleSrc  (&oracleMat0 ), * pOracleDest  (&oracleMat1);
  Matrix* pBidderSrc  (&bidderMat0 ), * pBidderDest  (&bidderMat1);
  // iteration vars
  auto search = make_search_engine();
  std::pair<double,double> maxPair, bestMeanInterval;
  std::pair<int, double> bidderKP, oracleKP;
  double oracleBid,bidderBid;
  bestMeanInterval = std::make_pair(10,0);
  int done = 1;
  std::cout << std::setprecision(3);  // debug
  for (int round = nRounds; round > 0; --round, ++done)
  { if (use0)
    { pUtilitySrc = &utilityMat0;    pOracleSrc  = &oracleMat0;   pBidderSrc  = &bidderMat0;
      pUtilityDest= &utilityMat1;    pOracleDest = &oracleMat1;   pBidderDest = &bidderMat1;
    } else
    { pUtilitySrc = &utilityMat1;    pOracleSrc  = &oracleMat1;   pBidderSrc  = &bidderMat1;
      pUtilityDest= &utilityMat0;    pOracleDest = &oracleMat0;   pBidderDest = &bidderMat0;
    }
    use0 = !use0;
    //    std::cout << " ---------  Round " << round << " ----------------------\n";
    // recursion is in 'reverse time' starting from final round, working forward
    // ko,kb denote positions in wealth array; higher values mean higher wealth
    // KP pairs indicate position and weight if reject
#ifdef  SHOW_PROGRESS
    const bool show_progress (true);
#else
    const bool show_progress (false);
#endif
    for (int ko=done; ko<nColumns-1; ++ko) 
    { oracleBid = oracleWealth.bid(ko);
      oracleKP  = oracleWealth.wealth_position(ko);   
      for (int kb=done; kb<nColumns-1; ++kb) 
      { bidderBid = bidderWealth.bid(kb);
	bidderKP  = bidderWealth.wealth_position(kb);
	utility.set_constants(oracleBid, bidderBid,                               // oracle is alpha, bidder is beta;  bidder position on cols
			      (*pUtilitySrc)(ko-1    , kb-1     ),                // v00  neither rejects
			      reject_value  (ko-1    , bidderKP, *pUtilitySrc),   // v01  only bidder rejects
			      reject_value  (oracleKP, kb-1    , *pUtilitySrc),   // v10  only oracle rejects
			      reject_value  (oracleKP, bidderKP, *pUtilitySrc));  // v11  both reject
	maxPair = search.find_maximum(utility);       // returns opt, f(opt)
	// monitor range of optimal means
	if(maxPair.first < bestMeanInterval.first)
	  bestMeanInterval.first = maxPair.first;
	else if (maxPair.first > bestMeanInterval.second)
	  bestMeanInterval.second = maxPair.first;
	double utilAtMuEqualZero = utility(0.0);
	if (maxPair.second < utilAtMuEqualZero)
	  maxPair = std::make_pair(0.0,utilAtMuEqualZero);
	// save mean if oracle in desired wealth state
	if      (mIndexA == ko) meanMatA(round+1,kb) = maxPair.first;
	else if (mIndexB == ko) meanMatB(round+1,kb) = maxPair.first; 
	(*pUtilityDest)(ko,kb) = maxPair.second;
	// oracle on rows of outcome, bidder on columns, fill nests down to lower right corner _|
	(* pOracleDest)(ko,kb) = utility.oracle_utility(maxPair.first,
							(*pOracleSrc)(ko-1    , kb-1   ),
							reject_value (ko-1    , bidderKP, *pOracleSrc, show_progress),
							reject_value (oracleKP, kb-1    , *pOracleSrc, show_progress),
							reject_value (oracleKP, bidderKP, *pOracleSrc, show_progress));  
        (* pBidderDest)(ko,kb) = utility.bidder_utility(maxPair.first,
							(*pBidderSrc)(ko-1    , kb-1   ),
							reject_value (ko-1    , bidderKP, *pBidderSrc),
							reject_value (oracleKP, kb-1    , *pBidderSrc),
							reject_value (oracleKP, bidderKP, *pBidderSrc));
	/* huge debugging output... */
#ifdef  SHOW_PROGRESS
	std::cout << std::endl
		  << " ****  [ ko=" << ko << ", kb=" << kb << "] "
		  << " oracle bid " << oracleBid << " with kp = (" << oracleKP.first << "," << oracleKP.second << ")   "
		  << " bidder bid " << bidderBid << " with kp = (" << bidderKP.first << "," << bidderKP.second << ")" << std::endl;
	std::cout << "   alpha=" << utility.alpha() << "  beta=" << utility.beta() 
		  << "  Oracle value is " << (*pOracleDest)(ko,kb) << "= Util(" << maxPair.first << " , "
		  << (*pOracleSrc)(ko-1    , kb-1   ) << " , "
		  << reject_value (ko-1    , bidderKP, *pOracleSrc) << " , "
		  << reject_value (oracleKP, kb-1    , *pOracleSrc) << " , "
		  << reject_value (oracleKP, bidderKP, *pOracleSrc) 
		  << ")   Bidder value is " << (*pBidderDest)(ko,kb) << "= Util(" << maxPair.first << " , "
		  << (*pBidderSrc)(ko-1    , kb-1   ) << " , "
		  << reject_value (ko-1    , bidderKP, *pBidderSrc) << " , "
		  << reject_value (oracleKP, kb-1    , *pBidderSrc) << " , "
		  << reject_value (oracleKP, bidderKP, *pBidderSrc) << ")" << std::endl;
#endif
      }
      // fill padding column
      (*pUtilityDest)(ko,nColumns-1) = (* pUtilityDest)(ko,nColumns-2);
      (*pOracleDest )(ko,nColumns-1) = (* pOracleDest )(ko,nColumns-2);
      (*pBidderDest )(ko,nColumns-1) = (* pBidderDest )(ko,nColumns-2);
    }
    // fill padding row
    for (int kb=0; kb<nColumns; ++kb) 
    { (*pUtilityDest)(nColumns-1,kb) = (* pUtilityDest)(nColumns-2,kb);
      (*pOracleDest )(nColumns-1,kb) = (* pOracleDest )(nColumns-2,kb);
      (*pBidderDest )(nColumns-1,kb) = (* pBidderDest )(nColumns-2,kb);
    }
#ifdef  SHOW_PROGRESS
    std::cout << "Based on using oracle source matrix " << std::endl;
    std::cout << *pOracleSrc << std::endl;
    std::cout << " - - - - - - - - - - - - - - - - -" << std::endl << std::endl;
#endif
  }
  std::cout << std::setprecision(6);
  if(writeDetails)
  { std::ostringstream ss;
    int angle (utility.angle());
    ss << "runs/bellman2.g" << angle << ".n" << nRounds << ".";
    write_matrix_to_file(ss.str() + "meanA"  ,    meanMatA  );
    write_matrix_to_file(ss.str() + "meanB"  ,    meanMatB  );
    bool writeOneRow = false;
    if (!writeOneRow) // write final destination matrices
    { write_matrix_to_file(ss.str() + "utility", *pUtilityDest);
      write_matrix_to_file(ss.str() + "oracle" ,  *pOracleDest);
      write_matrix_to_file(ss.str() + "bidder" ,  *pBidderDest);
    }
    else // write omega row of each to summary file
    { std::string fileName (ss.str() + "summary");
      write_vector_to_file(fileName, pUtilityDest->row(iOmega));
      write_vector_to_file(fileName,  pOracleDest->row(iOmega), true);  // append 
      write_vector_to_file(fileName,  pBidderDest->row(iOmega), true);
    }
  }
  // write summary of configuration and results to stdio
  std::cout << utility.angle() << " " << bidderWealth.omega() << "   " << nRounds   << "   "
	    << (*pUtilityDest)(nRounds,nRounds) << " " << (*pOracleDest)(nRounds,nRounds) << " " << (*pBidderDest)(nRounds,nRounds) << std::endl;
}





// --------------------------------------------------------------------------------------------------------------
//
//      Constrained    Constrained    Constrained    Constrained    Constrained    Constrained    Constrained    
//
//      This version handles state-dependent expert, returning the grid of expert and bidder values
//      needed for the next round of the recursion.
//
// --------------------------------------------------------------------------------------------------------------


void
solve_constrained_bellman_alpha_equation (double gamma, double omega, int nRounds, double spendPct, double oracleProb, Distribution const& bidderProb, bool writeDetails)
{
  const int maxIterations (200);   
  const double tolerance  (0.0001);
  const double grid       (0.5);
  const std::pair<double,double> searchInterval = std::make_pair(0.05,10.0);
  Line_Search::GoldenSection search(tolerance, searchInterval, grid, maxIterations);
  ConstrainedExpertCompetitiveAlphaGain compRatio (gamma, omega, spendPct, oracleProb, bidderProb);
    
  // space for holding intermediate results; extra row for boundary condition
  // code flips between these on read and write
  Matrix gain0   = Matrix::Constant (nRounds+1, nRounds+1, omega - gamma * omega);
  Matrix oracle0 = Matrix::Constant (nRounds+1, nRounds+1, omega);
  Matrix bidder0 = Matrix::Constant (nRounds+1, nRounds+1, omega);
  Matrix gain1   = Matrix::Zero (nRounds+1, nRounds+1);
  Matrix oracle1 = Matrix::Zero (nRounds+1, nRounds+1);
  Matrix bidder1 = Matrix::Zero (nRounds+1, nRounds+1);
  Matrix mean    = Matrix::Zero (nRounds, nRounds);
  bool use0 = true;

  // Check the probability distribution being used
  //  std::cout << std::endl << std::endl;
  //  for (int i=0; i<10; ++i) std::cout << i << " " << omega * spendPct * bidderProb(i,0) << "    ";
  //  std::cout << std::endl << std::endl;

  // alternate between reading and writing these matrices
  Matrix* pGainSrc;   Matrix* pGainDest = NULL;
  Matrix* pOracleSrc; Matrix* pOracleDest = NULL;
  Matrix* pBidderSrc; Matrix* pBidderDest = NULL;
  std::pair<double,double> bestMeanInterval(100,0);
  for (int round = nRounds; round > 0; --round)
  { if (use0)
    { pGainSrc    = &gain0;    pOracleSrc  = &oracle0;   pBidderSrc  = &bidder0;
      pGainDest   = &gain1;    pOracleDest = &oracle1;   pBidderDest = &bidder1;
    } else
    { pGainSrc    = &gain1;    pOracleSrc  = &oracle1;   pBidderSrc  = &bidder1;
      pGainDest   = &gain0;    pOracleDest = &oracle0;   pBidderDest = &bidder0;
    }
    use0 = !use0;
    if (writeDetails)
    { std::cout << "\n\n--------------- Round " << round << " comparative value source --------------------- " << std::endl << pGainSrc->topLeftCorner(round+1,round+1) << std::endl;
      std::cout << " --------------- Expert (top 2 rows) --------------------- " << std::endl << pOracleSrc->topLeftCorner(imin(2,round+1),round+1) << std::endl;
      std::cout << " --------------- Bidder (top 2 rows) --------------------- " << std::endl << pBidderSrc->topLeftCorner(imin(2,round+1),round+1) << std::endl;
    }
    for (int i=0; i<round; ++i)        // next round status of expert
    { for (int j=0; j<round; ++j)      //                      bidder
      { std::pair<double,double> maxPair;
	compRatio.set_delay (i, j, round, nRounds, (*pGainSrc)(0,0), (*pGainSrc)(i+1,0), (*pGainSrc)(0,j+1), (*pGainSrc)(i+1,j+1));
	maxPair = search.find_minimum(compRatio);        // mean, f(mean)
	if(maxPair.first < bestMeanInterval.first)       // monitor range of optimal means
	  bestMeanInterval.first = maxPair.first;
	else if (maxPair.first > bestMeanInterval.second)
	  bestMeanInterval.second = maxPair.first;
	double atZero = compRatio(0.0);
	if (maxPair.second < atZero)
	  maxPair = std::make_pair(0.0,atZero);
	mean(i,j) = maxPair.first;
	(*pGainDest)(i,j) = maxPair.second;
	(*pOracleDest)(i,j) = compRatio.value_to_oracle(maxPair.first, (*pOracleSrc)(0,0), (*pOracleSrc)(i+1,0), (*pOracleSrc)(0,j+1), (*pOracleSrc)(i+1,j+1));
	(*pBidderDest)(i,j) = compRatio.value_to_bidder(maxPair.first, (*pBidderSrc)(0,0), (*pBidderSrc)(i+1,0), (*pBidderSrc)(0,j+1), (*pBidderSrc)(i+1,j+1));
      }
    }
    if (writeDetails)
      std::cout << "\n\n---------------   MEAN  --------------------- " << std::endl << mean.topLeftCorner(round,round) << std::endl;
  }
  // write parms and final values to std io
  std::clog << "Interval for optimal means is [" << bestMeanInterval.first << "," << bestMeanInterval.second << std::endl;
  std::cout << gamma             << " " << omega               << " " << nRounds             << " " << spendPct << " "
	    << searchInterval.first << " " << searchInterval.second << " " 
	    << (*pGainDest)(0,0) << " " << (*pOracleDest)(0,0) << " " << (*pBidderDest)(0,0) << std::endl;
}

 /////////////////////////////////  Constrained Expert  ///////////////////////////////////////


 void
 ConstrainedExpertCompetitiveAlphaGain::set_delay (int i, int j, int /* t */, int /* nRounds */, double v00, double vi0, double v0j, double vij)
 {
   mAlpha = mOmega             * mExpertDist(i /*,nRounds-t*/ );  // no spending constraint for expert
   mBeta  = mOmega * mSpendPct * mBidderProb(j /*,nRounds-t*/ );
   mV00 = v00;
   mVi0 = vi0;
   mV0j = v0j;
   mVij = vij;
 }

 double
 ConstrainedExpertCompetitiveAlphaGain::operator()(double mu) const
 {
   double a = mAlpha;                 // a = optimal_alpha(mu, mOmega); 
   double ra = reject_prob(mu,a);
   double rb = reject_prob(mu,mBeta);
   double gain = (mOmega * ra - a) - mGamma * (mOmega * rb - mBeta) + ra * rb * mV00 + (1-ra) * (1-rb) * mVij + ra * (1-rb) * mV0j + (1-ra) * rb * mVi0;
   // mondo output
   // if (mu == 0.0) std::cout << "\n gain @ 0 = " << gain << " using  a,ra = " << a << "," << ra << "   b,rb = " << mBeta << "," << rb
   //      		      << "    v " << mV00 << " " << mVij << " " << mVi0 << " " << mV0j << std::endl;
   return gain;
 }


 double
 ConstrainedExpertCompetitiveAlphaGain::value_to_oracle(double mu, double v00, double vi0, double v0j, double vij) const
 {
   double value, ra, rb;
   if (mu < 0.00001)
   { ra = mAlpha;                          // 0 if unconstrained
     value = mAlpha * (mOmega-1.0);        // 0 if unconstrained
     rb = mBeta;
   }
   else
   { double a = mAlpha;                    // a=optimal_alpha(mu, mOmega);
     ra = reject_prob(mu,a);
     value = mOmega * ra - a ;
     rb = reject_prob(mu, mBeta);
   }
   return  value + ra * rb * v00 +  (1-ra) * rb * vi0 + ra * (1-rb) * v0j + (1-ra) * (1-rb) * vij;
 }


double
ConstrainedExpertCompetitiveAlphaGain::value_to_bidder(double mu, double v00, double vi0, double v0j, double vij) const
{
  double ra, rb;
  if (mu < 0.00001)
  { ra = mAlpha;                           //     would be zero if unconstrained
    rb = mBeta;
  }
  else
  { ra = reject_prob(mu,mAlpha);           //     ra = reject_prob(mu,optimal_alpha(mu, mOmega));
    rb = reject_prob(mu, mBeta);
  }
  return  (mOmega * rb - mBeta)  + ra * rb * v00 +  (1-ra) * rb * vi0 + ra * (1-rb) * v0j + (1-ra) * (1-rb) * vij;
}



// --------------------------------------------------------------------------------------------------------------
//
//      Unconstrained     Unconstrained     Unconstrained     Unconstrained     Unconstrained     Unconstrained 
//
// --------------------------------------------------------------------------------------------------------------


// this version does alpha-wealth rather than input utility
void
solve_bellman_alpha_equation (double gamma, double omega, int nRounds, double spendPct, Distribution const& f, bool writeDetails)
{
  const int maxIterations (100);
  const double tolerance  (0.0001);
  const double grid       (0.5);
  const std::pair<double,double> searchInterval = std::make_pair(1.5,6.5);
  Line_Search::GoldenSection search(tolerance, searchInterval, grid, maxIterations);
  ExpertCompetitiveAlphaGain compRatio (gamma, omega, f, spendPct);
     
  // space for holding intermediate results
  Matrix gain  (nRounds+1, nRounds+1);   // extra row for boundary condition
  Matrix oracle(nRounds+1, nRounds+1);
  Matrix bidder(nRounds+1, nRounds+1);
  Matrix mean = Matrix::Zero(nRounds,nRounds);
  // capture further results
  //  Matrix prob = Matrix::Zero(nRounds,nRounds);
  
  //  initialize boundary conditions
  double v0   = omega - gamma * omega;
  for (int j=0; j<nRounds+1; ++j)
  { gain(nRounds,j)  = v0;
    oracle(nRounds,j)=omega;
    bidder(nRounds,j)=omega;
  }
  // stores intermediates in rows of triangular array
  for (int row = nRounds-1; row > -1; --row)
  { v0 = gain(row+1,0);
    double b0 = bidder(row+1,0);
    double o0 = oracle(row+1,0);
    for (int col=0; col<=row; ++col)
    { std::pair<double,double> maxPair;
      compRatio.set_k(col, nRounds-1-row, v0, gain(row+1,col+1)); 
      double atZero = compRatio(0.0);
      maxPair = search.find_minimum(compRatio);
      if (maxPair.second < atZero)
	maxPair = std::make_pair(0.0,atZero);
      gain  (row,col) = maxPair.second;
      oracle(row,col) = compRatio.value_to_oracle(maxPair.first, o0, oracle(row+1,col+1));
      bidder(row,col) = compRatio.value_to_bidder(maxPair.first, b0, bidder(row+1,col+1));
      mean (row,col) = maxPair.first;
    }
  }
  // write solution (without boundary row) to file
  if(writeDetails)
  { std::string fileName;
    std::ostringstream ss;
    int gammaInt (trunc(10 * gamma));
    ss << "bellman.g" << gammaInt << ".n" << nRounds << ".";
    write_matrix_to_file(ss.str() + "gain",     gain.topLeftCorner(gain.rows()-1, gain.rows()-1));  // omit boundary row
    write_matrix_to_file(ss.str() + "oracle", oracle.topLeftCorner(oracle.rows()-1,oracle.rows()-1));
    write_matrix_to_file(ss.str() + "bidder", bidder.topLeftCorner(bidder.rows()-1, bidder.rows()-1));
    write_matrix_to_file(ss.str() + "mu",       mean.topLeftCorner(mean.rows(), mean.rows()));
  }
  // write parameters and final values to std io
  std::cout <<  gamma     << " " << omega       << " " << nRounds     << " " << spendPct << " "
	    << searchInterval.first << " " << searchInterval.second  << " " 
	    << gain(0,0) << " " << oracle(0,0) << " " << bidder(0,0) << std::endl;
}



/////////////////////////////////  Unconstrained Expert  ///////////////////////////////////////


double
ExpertCompetitiveAlphaGain::operator()(double mu) const
{
  if(mu < 0.00001)
  { // std::cout << "For mu = 0 reject prob r = b = " << b << " and V0 = " << mV0 << " and V_k+1 = " << mVkp1 << std::endl;
    return mGamma * mBetaK * (1.0-mOmega) + mBetaK*mV0 + (1-mBetaK)*mVkp1;
  }
  else
  {
    double rb = reject_prob(mu,mBetaK);
    double a = optimal_alpha(mu, mOmega);
    double ra = reject_prob(mu,a);
    double gain = (mOmega * ra - a) - mGamma * (mOmega * rb - mBetaK) + rb * mV0 + (1-rb) * mVkp1;
    return gain;
  }
}
  
double
ExpertCompetitiveAlphaGain::value_to_oracle(double mu, double o0, double okp1) const 
{
  double value, rb;
  if (mu < 0.00001)
  { value = 0.0;
    rb = mBetaK;
  }
  else
  { double a = optimal_alpha(mu, mOmega);
    double ra = reject_prob(mu,a);
    value = mOmega * ra - a ;
    rb = reject_prob(mu, mBetaK);
  }
  return value + rb * o0 + (1-rb) * okp1;
}
  
double
ExpertCompetitiveAlphaGain::value_to_bidder (double mu, double b0, double bkp1) const
{
  double rb = (mu < 0.00001) ? mBetaK : reject_prob(mu,mBetaK);
  return mOmega * rb - mBetaK + rb * b0 + (1-rb) * bkp1;
}

