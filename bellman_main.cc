#include "bellman.h"

#include "line_search.Template.h"
#include "wealth.Template.h"

#include <math.h>
#include <iostream>
#include <getopt.h>
#include "read_utils.h"     

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>


// Where to start the universal coder

const int universalStart (1);


// player probability and omega, i p o, (initial prob omega)
//    omega = 1 defines an unconstrained player (oracle); non-zero omega implies contrained
//    omega = 0 implies constant alpha level player (old school statistician)
//    prob = 1  risk inflation oracle
//    prob = 0  ls oracle

typedef boost::tuple<double, double, double> Triple;

double w0   (Triple const& p) { return boost::get<0>(p); }
double prob (Triple const& p) { return boost::get<1>(p); }
double omega(Triple const& p) { return boost::get<2>(p); }


// prob=0 signals universal, prob > 0 is geometric

DualWealthArray*
make_wealth_array(Triple const& parms, double scale, int nRounds);


// used to format output file names
int
round_parm(double x)
{
  return floor(100 * x);
}


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
  double     scale    = 1.0;                        // multiplier of universal code in unconstrained
  bool     writeTable = false;                      // if false, only return final value
  Triple    oracle    = boost::make_tuple(0,0,0);   //   (W0, alpha, oracle omega) set omega to nonzero to contrain; else unconstrained
  Triple    bidder    = boost::make_tuple(0,0,0);   //   (W0, beta, bidder omega)

  parse_arguments(argc, argv, riskUtil, angle, oracle, bidder, scale, nRounds, writeTable);

  /*
     Note that alpha (aka, the oracle probability for a Bayes oracle)
     'lives' in the utility function object, and W0 and omega are
     in the wealth function for encoding the index position of the
     wealth when a rejection occurs.
  */
  
  std::clog << "MAIN: Building bidder wealth array for " << nRounds << " rounds with (w0, beta, omega) =" << bidder << ", and scale=" << scale << std::endl;
  DualWealthArray *pBidderWealth = make_wealth_array(bidder, scale, nRounds);
  // pBidderWealth->write_to(std::clog, true); std::clog << std::endl; // as lines

  if(omega(oracle) == 1) // unconstrained competitor
  { std::cout << "Oracle(" << w0(oracle) << "," << prob(oracle) << "," << omega(oracle) << ") " << pBidderWealth->name() << " ";
    if (riskUtil)
    { RiskVectorUtility utility(angle, prob(oracle));
      solve_bellman_utility (nRounds, utility, *pBidderWealth, writeTable);
    }
    else
    { RejectVectorUtility utility(angle, prob(oracle));
      solve_bellman_utility (nRounds, utility, *pBidderWealth, writeTable);
    }
  }
  else                     // constrained competitor needs to track state as well
  { std::clog << "MAIN: Column player (bidder) has (w0,p,omega)=" << bidder << " and uses wealth array ... " << *pBidderWealth <<  std::endl;
    DualWealthArray *pOracleWealth = make_wealth_array(oracle, scale, nRounds);
    std::clog << "MAIN: Row player (oracle) has (w0,0,omega)=" << oracle << " and uses wealth array ... " << *pOracleWealth << std::endl;
    std::cout << pOracleWealth->name() << " "     << pBidderWealth->name() << " ";    // start of file output
    std::ostringstream ss;
    ss << angle << " " << scale << " "
       << prob(oracle) << " " << omega(oracle) << " "
       << prob(bidder) << " " << omega(bidder) << " ";
    if (riskUtil)
    { RiskMatrixUtility utility(angle);
      ss << ".risk";
      solve_bellman_utility (nRounds, utility, *pOracleWealth, *pBidderWealth, ss.str(), writeTable);
    }
    else
    { RejectMatrixUtility utility(angle);
      ss << ".reject";
      solve_bellman_utility (nRounds, utility, *pOracleWealth, *pBidderWealth, ss.str(), writeTable);
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
  bool setBidderW0 = false;
  bool setOracleW0 = false;
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
	setOracleW0 = true;
	boost::get<0>(oracleIPO) = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'o' : 
      {
	boost::get<1>(oracleIPO) = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'O' : 
      {
	boost::get<2>(oracleIPO) = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'I' : 
      {
	setBidderW0 = true;
	boost::get<0>(bidderIPO) = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'b' : 
      {
	boost::get<1>(bidderIPO) = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'B' : 
      {
	boost::get<2>(bidderIPO) = read_utils::lexical_cast<double>(optarg);
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
  if(!setBidderW0)  boost::get<0>(bidderIPO) = omega(bidderIPO);
  if(!setOracleW0)  boost::get<0>(oracleIPO) = omega(oracleIPO);
}


DualWealthArray*
make_wealth_array(Triple const& parms, double scale, int nRounds)
{
  if (0 == omega(parms))             // fixed alpha bidder
  { double p = prob(parms)/nRounds;
    std::clog << "MAIN: Fixed alpha bidder with constant bid alpha=" << p << std::endl;
    return new DualWealthArray(p);
  }
  else if(0 == prob(parms))          // universal
  { std::clog << "MAIN: Making universal array with scale=" << scale << " and W0=" << w0(parms) << " omega=" << omega(parms) << std::endl;
    UniversalBidder bidder(scale);
    return new DualWealthArray(bidder.identifier(), w0(parms), omega(parms), bidder, nRounds);
  }
  else 
  { std::clog << "MAIN: Making geometric wealth array with p=" << prob(parms) << ", scale=" << scale
	      << " and W0=" << w0(parms) << " omega=" << omega(parms) << std::endl;
    UniversalBidder univ(scale);
    GeometricBidder geoBidder(prob(parms), univ.total_wealth());
    return new DualWealthArray(geoBidder.identifier(), w0(parms), omega(parms), geoBidder, nRounds);
  }
}
