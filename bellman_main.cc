#include "bellman.h"

#include "line_search.Template.h"
#include "wealth.Template.h"

#include <math.h>
#include <iostream>
#include <getopt.h>
#include "read_utils.h"     


// Where to start the universal coder

const int universalStart (1);

// Need to use special order of calls to fill geometric
// prob=0 signals universal, prob > 0 is geometric

DualWealthArray*
make_wealth_array(int nRounds, double omega, double prob, double scale);

int
round_parm(double x)
{
  return floor(100 * x);
}

//  prob character indicates the distribution, u for universal and g for geometric

void
parse_arguments(int argc, char** argv,
		bool &riskUtil, double &angle, int &nRounds, bool &constrained,
		double &oracleProb, double &bidderProb,  double &scale,
		double &omega, bool &writeTable);


// Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     
int  main(int argc, char** argv)
{
  // default arguments
  bool      riskUtil  = false;    // risk or rejection, default is rejection (which is fast)
  double       angle  =     0;    // in degrees
  int        nRounds  =   100;
  bool     constrain  = false;    // ignores oracle prob if not constrained
  double   oracleProb = 0.0;      // use univ oracle if zero prob and constrained
  double   bidderProb = 0.0;
  double     scale    = 1.0;      // multiplier of universal code in unconstrained
  bool     writeTable = false;    // if false, only return final value
  double     omega    = 0.05;     // also sets the initial wealth

  parse_arguments(argc, argv, riskUtil, angle, nRounds, constrain, oracleProb, bidderProb, scale, omega, writeTable);

  /*
     Note that alpha (aka, the oracle probability for a Bayes oracle)
     'lives' in the utility function object, and omega is embedded
     into the wealth function for encoding the index position of the
     wealth when a rejection occurs.
  */
  
  std::clog << "MAIN: Building wealth array for " << nRounds << " rounds with omega=" << omega
	    << ", bidder prob=" << bidderProb << ", and scale=" << scale << std::endl;
  DualWealthArray *pBidderWealth = make_wealth_array(nRounds, omega, bidderProb, scale);
  // pBidderWealth->write_to(std::clog, true); std::clog << std::endl; // as lines
  
  if(!constrain)           // unconstrained competitor
  { std::cout << "Oracle(" << oracleProb << ") " << pBidderWealth->name() << " ";
    if (riskUtil)
    { RiskVectorUtility utility(angle, oracleProb);
      solve_bellman_utility (nRounds, utility, *pBidderWealth, writeTable);
    }
    else
    { RejectVectorUtility utility(angle, oracleProb);
      solve_bellman_utility (nRounds, utility, *pBidderWealth, writeTable);
    }
  }
  else                     // constrained competitor needs to track state as well
  { std::clog << "MAIN: Column player (bidder) uses... " << *pBidderWealth << std::endl;
    DualWealthArray *pOracleWealth = make_wealth_array(nRounds, omega, oracleProb, scale);
    std::clog << "MAIN: Row player (oracle) uses wealth array " << *pOracleWealth << std::endl;
    std::cout << pOracleWealth->name() << " "     << pBidderWealth->name() << " ";
    std::ostringstream ss;
    ss << "druns/bellman.a" << angle << ".s" << round_parm(scale) <<".o" << round_parm(omega) << ".op" << round_parm(oracleProb) << ".bp" << round_parm(bidderProb);
    if (riskUtil)
    { RiskMatrixUtility utility(angle, omega);
      ss << ".risk";
      solve_bellman_utility (nRounds, utility, *pOracleWealth, *pBidderWealth, ss.str(), writeTable);
    }
    else
    { RejectMatrixUtility utility(angle, omega);
      ss << ".reject";
      solve_bellman_utility (nRounds, utility, *pOracleWealth, *pBidderWealth, ss.str(), writeTable);
    }
  }
  return 0;
}



void
parse_arguments(int argc, char** argv,
		bool &riskUtil, double &angle, int &nRounds, bool &constrain,
		double &oracleProb, double &bidderProb,   // zero denotes universal
		double &scale, double &omega, bool &writeTable)
{
  static struct option long_options[] = {
    {"risk",             no_argument, 0, 'R'},
    {"reject",           no_argument, 0, 'r'},
    {"angle",      required_argument, 0, 'a'},
    {"constrain",        no_argument, 0, 'c'},
    {"oracleprob", required_argument, 0, 'o'},
    {"bidderprob", required_argument, 0, 'b'},
    {"scale",      required_argument, 0, 's'},
    {"rounds",     required_argument, 0, 'n'},
    {"omega",      required_argument, 0, 'W'},
    {"write",            no_argument, 0, 'w'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  bool rejectUtil = true;
  while (-1 !=(key = getopt_long (argc, argv, "Rra:co:b:s:n:W:w", long_options, &option_index))) // colon means has argument
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
    case 'c' : 
      {
	constrain = true;
	break;
      }
    case 'o' : 
      {
	oracleProb = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'b' : 
      {
	bidderProb = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 's' : 
      {
	scale = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'W' :
      {
	omega = read_utils::lexical_cast<double>(optarg);
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
}


DualWealthArray*
make_wealth_array(int nRounds, double omega, double prob, double scale)
{
  if (0 == omega)             // fixed alpha bidder
  { prob = prob/nRounds;
    std::clog << "MAIN: Fixed alpha bidder with constant bid alpha=" << prob << std::endl;
    return new DualWealthArray(prob);
  }
  else if(0 == prob)          // universal
  { std::clog << "MAIN: Making universal array with scale=" << scale << " and omega=" << omega << std::endl;
    UniversalBidder bidder(scale);
    return new DualWealthArray(bidder.identifier(), omega, omega, bidder, nRounds);
    //  return new WealthArray(omega, omega, nRounds, ScaledUniversalDist(scale));
  }
  else 
  { std::clog << "MAIN: Making geometric array with prob=" << prob << ", scale=" << scale << " and omega=" << omega << std::endl;
    UniversalBidder univ(scale);
    GeometricBidder geoBidder(prob, univ.total_wealth());
    return new DualWealthArray(geoBidder.identifier(), omega, omega, geoBidder, nRounds);
  }
}
