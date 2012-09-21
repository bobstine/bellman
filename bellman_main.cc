#include "bellman.h"

#include <math.h>

#include <iostream>
#include <sstream>
#include <getopt.h>
#include "read_utils.h"     


// Where to start the universal coder

const int universalStart (1);

// Need to use special order of calls to fill geometric
// prob=0 signals universal, prob > 0 is geometric

WealthArray*
make_wealth_array(double omega, int iOmega, double prob);


//  prob character indicates the distribution, u for universal and g for geometric

void
parse_arguments(int argc, char** argv,
		bool &riskUtil, double &angle, int &nRounds, bool &constrained,
		double &oracleProb, double &bidderProb,  
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
  bool     writeTable = false;    // if false, only return final value
  double     omega    = 0.05;     // also sets the initial wealth

  parse_arguments(argc, argv, riskUtil, angle, nRounds, constrain, oracleProb, bidderProb, omega, writeTable);
  const int iOmega    (nRounds+1);
  
  WealthArray* pBidderWealth = make_wealth_array(omega, iOmega, bidderProb);

  if(!constrain)           // unconstrained oracle 
  { std::cout << "uncon(" << oracleProb << ") " << pBidderWealth->name() << " ";
    if (riskUtil)
    { RiskVectorUtility utility(angle, omega);
      solve_bellman_utility (nRounds, utility, *pBidderWealth, writeTable);
    }
    else
    { RejectVectorUtility utility(angle, omega);
      solve_bellman_utility (nRounds, utility, *pBidderWealth, writeTable);
    }
  }
  else                     // constrained oracle needs wealth to track
  { WealthArray* pOracleWealth = make_wealth_array(omega, iOmega, oracleProb);
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
		double &omega, bool &writeTable)
{
  static struct option long_options[] = {
    {"risk",             no_argument, 0, 'R'},
    {"reject",           no_argument, 0, 'r'},
    {"angle",      required_argument, 0, 'a'},
    {"constrain",        no_argument, 0, 'c'},
    {"oracleprob", required_argument, 0, 'o'},
    {"bidderprob", required_argument, 0, 'b'},
    {"rounds",     required_argument, 0, 'n'},
    {"omega",      required_argument, 0, 'W'},
    {"write",            no_argument, 0, 'w'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  bool rejectUtil = true;
  while (-1 !=(key = getopt_long (argc, argv, "Rra:co:b:n:W:w", long_options, &option_index))) // colon means has argument
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
make_wealth_array(double omega, int iOmega, double prob)
{
  if(0 == prob)         // universal
  { double scale (2.0);
    std::clog << "MAIN: Making high-wealth universal array" << std::endl;
    return new WealthArray(omega, omega, iOmega, ScaledUniversalDist(scale));
  }
    // return new WealthArray(omega, iOmega, UniversalDist(universalStart));
  else if (prob > 1)    // uniform
    return new WealthArray(omega, iOmega, UniformDist( trunc(prob) ));
  else                  // geometric
    return new WealthArray(omega, iOmega, prob);
}
