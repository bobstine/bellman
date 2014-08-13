#include "bellman.Template.h"
#include "line_search.Template.h"
#include "wealth.Template.h"
#include "utility.Template.h"

#include <math.h>
#include <tuple>
#include <iostream>
#include <getopt.h>
#include "read_utils.h"     


// Where to start the universal coder

const int universalStart (1);

// player probability and omega, i p o, (initial prob omega)
//    omega = 1 defines an unconstrained player (oracle); non-zero omega implies contrained
//    omega = 0 implies constant alpha level player (old school statistician)
//    prob = 1  risk inflation oracle
//    prob = 0  ls oracle

typedef std::tuple<double, double, double> Triple;


double W0   (Triple const& p) { return std::get<0>(p); }
double prob (Triple const& p) { return std::get<1>(p); }
double omega(Triple const& p) { return std::get<2>(p); }

std::ostream & operator<<(std::ostream & o, Triple const& tri)
{ o << " { W0=" << W0(tri) << ", p=" << prob(tri) << ", w=" << omega(tri) << " } ";
  return o;
}


// prob=0 signals universal, prob > 0 is geometric
DualWealthArray*
make_wealth_array(Triple const& parms, double scale, int nRounds);


// used to format output file names
int
round_parm(double x)
{
  return floor(100 * x);
}

// read command line
void
parse_arguments(int argc, char** argv,
		bool &riskUtil, double &angle, Triple &oracle, Triple &bidder,
		double &scale, int &nRounds,  bool &writeTable);



// Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     
int  main(int argc, char** argv)
{
  // default arguments
  bool      riskUtil  = false;    // risk or rejection, default is rejection (which is fast)
  double       angle  =     0;    // in degrees
  int        nRounds  =   100;
  double     scale    = 1.0;                           // multiplier of universal code in unconstrained  (no longer used)
  bool     writeTable = false;                         // if false, only return final value
  Triple    oracle    = std::make_tuple(-1,-1,-1);   //   (W0, alpha, oracle omega) omega=1 implies unconstrained
  Triple    bidder    = std::make_tuple(-1,-1,-1);   //   (W0, beta, bidder omega)  negative values on exit parse were not set

  parse_arguments(argc, argv, riskUtil, angle, oracle, bidder, scale, nRounds, writeTable);

  std::clog << "MAIN: Running " << nRounds << " rounds at angle " << angle << " with writeTable=" << writeTable << std::endl;
  /*
     Note that alpha (aka, the oracle probability for a Bayes oracle)
     'lives' in the utility function object, and W0 and omega are
     in the wealth function for encoding the index position of the
     wealth when a rejection occurs.
  */
  
  std::clog << "MAIN: Building bidder wealth array for "
	    << nRounds << " rounds with " << bidder << ", and scale=" << scale << std::endl;
  DualWealthArray *pBidderWealth = make_wealth_array(bidder, scale, nRounds);
  // pBidderWealth->write_to(std::clog, true); std::clog << std::endl; // as lines

  if(omega(oracle) == 1)  // unconstrained competitor
  { std::clog << "MAIN: Oracle(W0,p,w)=" << oracle << " with bidder " << bidder << " and wealth function " << pBidderWealth->name() << std::endl;
    if (riskUtil)
    { RiskVectorUtility utility(angle, prob(oracle));
      solve_bellman_vector_utility (nRounds, utility, *pBidderWealth, writeTable);
    }
    else
    { RejectVectorUtility utility(angle, prob(oracle));
      solve_bellman_vector_utility (nRounds, utility, *pBidderWealth, writeTable);
    }
  }
  else                    // constrained competitor needs to track state as well
  { std::clog << "MAIN: Column player (bidder) " << bidder << " with wealth array ... " << *pBidderWealth <<  std::endl;
    DualWealthArray *pOracleWealth = make_wealth_array(oracle, scale, nRounds);
    std::clog << "MAIN: Row player (oracle)    " << oracle << " with wealth array ... " << *pOracleWealth << std::endl;
    std::clog << "MAIN: Players are : " << pOracleWealth->name() << " and " << pBidderWealth->name() << std::endl;
    std::ostringstream ss;
    ss << "n_" << nRounds << "_angle_" << angle << "_oracle_" << prob(oracle) << "_" << omega(oracle) << "_bidder_" << prob(bidder) << "_" << omega(bidder);
    AngleCriterion ac(angle);
    if (riskUtil)
    { RiskMatrixUtility<AngleCriterion> utility(ac);
      if (writeTable)
	solve_bellman_matrix_utility (nRounds, utility, *pOracleWealth, *pBidderWealth, ss.str(), writeTable); //tensor version
      else
	solve_bellman_matrix_utility (nRounds, utility, *pOracleWealth, *pBidderWealth);	
    }
    else
    { RejectMatrixUtility<AngleCriterion> utility(ac);
      solve_bellman_matrix_utility (nRounds, utility, *pOracleWealth, *pBidderWealth, ss.str(), writeTable);
    }
  }
  return 0;
}



void
parse_arguments(int argc, char** argv,
		bool &riskUtil, double &angle, Triple &oracleIPO, Triple &bidderIPO,
		double &scale, int &nRounds,  bool &writeTable)
{
  static struct option long_options[] = {
    {"risk",               no_argument, 0, 'R'},
    {"reject",             no_argument, 0, 'r'},
    {"angle",        required_argument, 0, 'a'},
    {"oracle_w0",    required_argument, 0, 'i'},
    {"oracle_prob",  required_argument, 0, 'o'},
    {"oracle_omega", required_argument, 0, 'O'},
    {"bidder_w0",    required_argument, 0, 'I'},
    {"bidder_prob",  required_argument, 0, 'b'},
    {"bidder_omega", required_argument, 0, 'B'},
    {"scale",        required_argument, 0, 's'},
    {"rounds",       required_argument, 0, 'n'},
    {"write",              no_argument, 0, 'w'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  bool rejectUtil = true;
  while (-1 !=(key = getopt_long (argc, argv, "Rra:i:o:O:I:b:B:s:n:w", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'R' : 
      {
	riskUtil = true; rejectUtil = false;
	break;
      }
    case 'r' : 
      {
	rejectUtil = true; riskUtil = false;
	break;
      }
    case 'a' : 
      {
	angle = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'n' :
      {
	nRounds = read_utils::lexical_cast<int>(optarg);
	break;
      }
    case 'i' : 
      {
	std::get<0>(oracleIPO) = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'o' : 
      {
	std::get<1>(oracleIPO) = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'O' : 
      {
	std::get<2>(oracleIPO) = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'I' : 
      {
	std::get<0>(bidderIPO) = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'b' : 
      {
	std::get<1>(bidderIPO) = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'B' : 
      {
	std::get<2>(bidderIPO) = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 's' : 
      {
	scale = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'w' : 
      {
	writeTable=true ;
	break;
      }
    default:
      {
	std::cout << "PARSE: Option not recognized; returning.\n";
      }
    } // switch
  } // while
  riskUtil = !rejectUtil;
  // set W0 to omega unless other value supplied
  if(W0(bidderIPO) < 0) std::get<0>(bidderIPO) = omega(bidderIPO);
  if(W0(oracleIPO) < 0) std::get<0>(oracleIPO) = omega(oracleIPO);
}


DualWealthArray*
make_wealth_array(Triple const& parms, double scale, int nRounds)
{
  const double maxWealth (5.0);
  
  if (1 == omega(parms))             // unconstrained, fixed wealth testimator
  { double w0 = W0(parms);
    std::clog << "MAIN: Fixed bidder with constant wealth=" << w0 << std::endl;
    return new DualWealthArray(w0);
  }
  else if(0 == prob(parms))          // universal
  { std::clog << "MAIN: Making universal array with scale=" << scale << " and W0=" << W0(parms) << " omega=" << omega(parms) << std::endl;
    UniversalRule bidder;
    return new DualWealthArray(bidder.identifier(), maxWealth, W0(parms), omega(parms), bidder, nRounds);
  }
  else                               // geometric
  { std::clog << "MAIN: Making geometric wealth array with p=" << prob(parms) << ", scale=" << scale
	      << " and W0=" << W0(parms) << " omega=" << omega(parms) << std::endl;
    GeometricRule geoBidder(prob(parms));
    return new DualWealthArray(geoBidder.identifier(), maxWealth, W0(parms), omega(parms), geoBidder, nRounds);
  }
}
