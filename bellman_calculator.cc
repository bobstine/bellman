/*
  Use this main routine to generate the two-point paths within the
  feasible set that explore the use of a Bayesian model to generate
  the feasible set.
*/

#include "bellman.h"
#include "wealth.Template.h"
#include "read_utils.h"  // coerce

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <ios>


inline std::ostream&
operator<<(std::ostream &output,  std::pair<double,double> const& pair)
{ output << " < " << pair.first << " , " << pair.second << " > "; return output; }




void
parse_arguments(int argc, char** argv,
		double &alpha, double &beta, double &omega, double &scale, int &nRounds, 
		double &signal);


// Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     
int  main(int argc, char** argv)
{
  const double maxWealth (5.0);

  // set default control values
  int    nRounds (100   );
  double alpha   (  0.1 );
  double beta    (  0.1 );
  double omega   (  0.5 );
  double scale   (  2.0 );
  double signal  (  2.5 );
  parse_arguments(argc, argv, alpha, beta, omega, scale, nRounds, signal);

  std::clog << "MAIN: settings are n=" << nRounds << " scale=" << scale << " omega=" << omega << " mu=" << signal << std::endl;
    
  // set up utility and wealth; angle 0.0 is not used
  UniversalRule bidder;
  DualWealthArray bidderWealth("Univ", maxWealth, omega, omega, bidder, nRounds);
  std::clog << "MAIN: Wealth array... " << bidderWealth << std::endl;
  RiskVectorUtility utility(0.0, alpha);

  // open output file if want further output, otherwise direct to stdio
  if (false) {
    std::ostringstream ss;
    ss << "riskpath_alpha" << alpha << "_beta" << beta << "_omega" << omega << "_scale" << scale << "_n" << nRounds<< "_mu" << signal ;
    std::ios_base::openmode mode = std::ios_base::trunc;
    std::ofstream output (ss.str(), mode);
  }
  // find expected risk, varying p_zero from very small to larger; write to output file
  { std::pair<double,double> result;
    double p;   
    p = 0.000001;
    result = find_process_risk (nRounds, p, signal, utility, bidderWealth);
    std::cout << p << " " << result.first << " " << result.second << std::endl;
    result = find_process_risk (nRounds, 1-p, signal, utility, bidderWealth);
    std::cout << 1-p << " " << result.first << " " << result.second << std::endl;
    p = 0.00001;
    result = find_process_risk (nRounds, p, signal, utility, bidderWealth);
    std::cout << p << " " << result.first << " " << result.second << std::endl;
    result = find_process_risk (nRounds, 1-p, signal, utility, bidderWealth);
    std::cout << 1-p << " " << result.first << " " << result.second << std::endl;
    p = 0.0001;
    result = find_process_risk (nRounds, p, signal, utility, bidderWealth);
    std::cout << p << " " << result.first << " " << result.second << std::endl;
    result = find_process_risk (nRounds, 1-p, signal, utility, bidderWealth);
    std::cout << 1-p << " " << result.first << " " << result.second << std::endl;
  }
  for (double p=0.001; p < .01; p+=0.001)
  { std::pair<double,double> result (find_process_risk (nRounds, p, signal, utility, bidderWealth));
    std::cout << p << " " << result.first << " " << result.second << std::endl;
    result = find_process_risk (nRounds, 1.0-p, signal, utility, bidderWealth);
    std::cout << 1.0-p << " " << result.first << " " << result.second << std::endl;
  }
  for (double p=0.01; p < 1; p+=0.01)
  { std::pair<double,double> result (find_process_risk (nRounds, p, signal, utility, bidderWealth));
    std::cout << p << " " << result.first << " " << result.second << std::endl;
  }
  
  return 0;
}


///////////////////

void
parse_arguments(int argc, char** argv,
		double &alpha,  double &beta, double &omega, double &scale, int &nRounds,
		double &signal)
{
  static struct option long_options[] = {
    {"alpha",      required_argument, 0, 'a'},
    {"beta" ,      required_argument, 0, 'b'},
    {"scale",      required_argument, 0, 's'},
    {"rounds",     required_argument, 0, 'n'},
    {"signal",     required_argument, 0, 'S'},
    {"omega",      required_argument, 0, 'W'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "a:b:s:n:N:p:S:W:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'a' : 
      {
	alpha = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'b' : 
      {
	beta = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 's' : 
      {
	scale = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'n' :
      {
	nRounds = read_utils::lexical_cast<int>(optarg);
	break;
      }
    case 'S' : 
      {
	signal = read_utils::lexical_cast<double>(optarg);
	break;
      }
    case 'W' :
      {
	omega = read_utils::lexical_cast<double>(optarg);
	break;
      }
    default:
      {
	std::cout << "PARSE: Option not recognized; returning.\n";
      }
    } // switch
  } // while
}

