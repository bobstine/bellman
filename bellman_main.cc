#include "bellman.h"

#include <math.h>

#include <iostream>
#include <getopt.h>
#include "read_utils.h"     


// Where to start the universal coder

const int universalStart (1);

// Need to use special order of calls to fill geometric
// prob=0 signals universal, prob > 0 is geometric

WealthArray*
make_wealth_array(int nRounds, double omega, double prob, double scale);


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
     'lives' in the utility function object, and omega is embedding
     into the wealth function for encoding the index position of the
     wealth when a rejection occurs.
  */
  
  std::clog << "MAIN: Building wealth array for " << nRounds << " rounds with omega=" << omega
	    << ", bidder prob=" << bidderProb << ", and scale=" << scale << std::endl;

  // WealthArray *pBidderWealth = make_wealth_array(nRounds, omega, bidderProb, scale);
  UniversalBidder bidder(scale);
  DualWealthArray *pBidderWealth = new DualWealthArray(bidder.identifier(), omega, omega, bidder, nRounds);
  std::clog << "MAIN: Bidder wealth array... " << *pBidderWealth << std::endl;
  
  if(!constrain)           // unconstrained oracle 
  { std::cout << "uncon(" << oracleProb << ") " << pBidderWealth->name() << " ";
    if (riskUtil)
    { RiskVectorUtility utility(angle, oracleProb);
      solve_bellman_utility (nRounds, utility, *pBidderWealth, writeTable);
    }
    else
    { RejectVectorUtility utility(angle, oracleProb);
      solve_bellman_utility (nRounds, utility, *pBidderWealth, writeTable);
    }
  }
  else                     // constrained oracle needs wealth to track
  { WealthArray* pOracleWealth = make_wealth_array(nRounds, omega, oracleProb, scale);
    std::cout << pOracleWealth->name() << " "     << pBidderWealth->name() << " ";
    if (riskUtil)
    { RiskMatrixUtility utility(angle, omega);
      solve_bellman_utility (nRounds, utility, *pOracleWealth, *pBidderWealth, writeTable);
    }
    else
    { RejectMatrixUtility utility(angle, omega); 
      solve_bellman_utility (nRounds, utility, *pOracleWealth, *pBidderWealth, writeTable);
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


WealthArray*
make_wealth_array(int nRounds, double omega, double prob, double scale)
{
  if(0 == prob)         // universal
  { std::clog << "MAIN: Making high-wealth universal array with scale=" << scale << " and omega=" << omega << std::endl;
    return new WealthArray(omega, omega, nRounds, ScaledUniversalDist(scale));
  }
  // return new WealthArray(omega, iOmega,UniversalDist(universalStart));
  else
  { int iZero (15);  // padding above initial wealth at omega
    if (prob > 1)    // uniform
      return new WealthArray(omega, iZero, nRounds, UniformDist( trunc(prob) ));
    else                  // geometric
      return new WealthArray(omega, iZero, nRounds, prob);
  }
}
