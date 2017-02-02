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
{ o << "{ W0=" << W0(tri) << ", p=" << prob(tri) << ", w=" << omega(tri) << "}";
  return o;
}


// prob=0 signals universal, prob > 0 is geometric
DualWealthArray*
make_wealth_array(Triple const& parms,  int nRounds);


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
		 int &nRounds,  bool &writeTable);



// Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     
int  main(int argc, char** argv)
{
  // default arguments
  double  angle = 0;
  
  bool      riskUtil  = false;    // risk or rejection, default is rejection (which is fast)
  //  double       RiB0   =  4.0;     // risk inflation intercept and slope
  double       RiB1   =  2.0;
  int        nRounds  = 200 ;
  bool     writeTable = false;                             // if false, only return final value
  Triple    oracle    = std::make_tuple(0.25, 0, 0.25);    //   (W0, univ, omega)              omega=1 implies unconstrained
  Triple    baseBidder= std::make_tuple(0.25,-1, 0.25);    //   (W0, beta, bidder omega)       negative values on exit parse were not set

  parse_arguments(argc, argv, riskUtil, angle, oracle, baseBidder,  nRounds, writeTable);

  std::clog << "MAIN: Oracle  " << oracle << std::endl;
  DualWealthArray *pOracleWealth = make_wealth_array(oracle,  nRounds);

  std::vector<double> psiVec = {.0001, 0.001, 0.01, 0.05, 0.10, 0.20, 0.30, 0.50};

  for(auto psi : psiVec)
  { Triple bidder = std::make_tuple(W0(baseBidder), psi, omega(baseBidder));
    std::cout << "MAIN: Bidder " << bidder << std::endl;
    DualWealthArray *pBidderWealth = make_wealth_array(bidder,  nRounds);
    RiskInflationCriterion ri(RiB1);
    RiskMatrixUtility<RiskInflationCriterion> utility(ri);
    solve_bellman_matrix_utility (nRounds, utility, *pOracleWealth, *pBidderWealth, " ", writeTable);
  }
  return 0;
}



void
parse_arguments(int argc, char** argv,
		bool &riskUtil, double &angle, Triple &oracleIPO, Triple &bidderIPO,
		 int &nRounds,  bool &writeTable)
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
    {"rounds",       required_argument, 0, 'n'},
    {"write",              no_argument, 0, 'w'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  bool rejectUtil = true;
  while (-1 !=(key = getopt_long (argc, argv, "Rra:i:o:O:I:b:B:n:w", long_options, &option_index))) // colon means has argument
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
make_wealth_array(Triple const& parms,  int nRounds)
{
  const double maxWealth (5.0);
  
  if (1 == omega(parms))             // unconstrained, fixed wealth testimator
  { double w0 = W0(parms);
    std::clog << "MAIN: Fixed bidder with constant wealth=" << w0 << std::endl;
    return new DualWealthArray(w0);
  }
  else if(0 == prob(parms))          // universal
  { std::clog << "MAIN: Making universal array with " << " W0=" << W0(parms) << " omega=" << omega(parms) << std::endl;
    UniversalRule bidder;
    return new DualWealthArray(bidder.identifier(), maxWealth, W0(parms), omega(parms), bidder, nRounds);
  }
  else                               // geometric
  { std::clog << "MAIN: Making geometric wealth array with p=" << prob(parms) 
	      << " and W0=" << W0(parms) << " omega=" << omega(parms) << std::endl;
    GeometricRule geoBidder(prob(parms));
    return new DualWealthArray(geoBidder.identifier(), maxWealth, W0(parms), omega(parms), geoBidder, nRounds);
  }
}
