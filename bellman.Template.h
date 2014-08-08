#include "bellman.h"
#include "utility.h"
#include "line_search.h"

#include <iomanip>
#include <ios>

#include "eigen_utils.h"
using EigenUtils::write_matrix_to_file;

//
//  Dual randomized wealth version
//

Line_Search::GoldenSection make_search_engine(void);

template<class Util>
void
solve_bellman_matrix_utility (int nRounds, Util &utility,
			      DualWealthArray const& rowWealth,  DualWealthArray const& colWealth, std::string config, bool writeTable)
{
  const int nRows (1+rowWealth.number_wealth_positions());                // extra 1 for padding; allow 0 * 0
  const int nCols (1+colWealth.number_wealth_positions());
  const std::pair<int,int> zeroIndex(std::make_pair(rowWealth.zero_index(), colWealth.zero_index()));
  std::clog << "BELL: Zero indices are " << zeroIndex.first << " & " << zeroIndex.second << std::endl;
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
	    << utility.identifier() << " " << 1000*rowWealth.omega()+colWealth.omega() << "   " << nRounds   << "   " 
	    << (*pUtilityDest)(zeroIndex.first, zeroIndex.second) << " "
	    << (*pRowDest    )(zeroIndex.first, zeroIndex.second) << " "
	    << (*pColDest    )(zeroIndex.first, zeroIndex.second) << std::endl;
}
