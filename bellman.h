#include <assert.h>

#include "utility.h"

#include <iostream>      // debug

/***********************************************************************************

  Bellman backward recursion for comptitive alpha-investing using
  several possible probability distributions, with spending percentage.

***********************************************************************************/

//  Finds the expected risk for process with probability p_0 for 0 and 1-p_0 for the given mean

std::pair<double,double>
  find_process_risk (int nRounds, double pZero, double mu, VectorUtility & utility, DualWealthArray const& bidderWealth);


//  These use a discrete wealth array to track the wealth of the bidder and (in constrained case) the oracle.
//  Both use a convex mixture of states when new wealth is not element of the array
//  Note: It's evil to pass in the reference,  but we don't care that the utility is modifiable; its there to be used.


// oracle with no wealth constraint
void
solve_bellman_vector_utility  (int nRounds, VectorUtility &util,                               DualWealthArray const& wealth, bool writeDetails);


// constrained oracle, two-player competition  (empty prefix means don't write)

//  this version does not save path information (alternating arrays)
template<class MatrixUtil>
void
solve_bellman_matrix_utility  (int nRounds, MatrixUtil & util,
			DualWealthArray const& oWealth,  DualWealthArray const& bidderWealth);

//  this version does save path information and writes out if asked (saves tensor)
template<class MatrixUtil>
void
solve_bellman_matrix_utility  (int nRounds, MatrixUtil & util,
			DualWealthArray const& oWealth,  DualWealthArray const& bidderWealth, std::string fileId, bool write);

