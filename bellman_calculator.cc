#include "bellman.h"
#include "read_utils.h"  // coerce

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <ios>


inline std::ostream&
operator<<(std::ostream &output, const std::pair<double,double> pair)
{ output << " < " << pair.first << " , " << pair.second << " > "; return output; }




void
parse_arguments(int argc, char** argv,
		int &nRounds, double &alpha,
		double &scale, double &omega, double &signal);


// Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     
int  main(int argc, char** argv)
{
  // set default control values
  int    nRounds (100   );
  double omega   (  0.5 );
  double alpha   (  0.05);
  double scale   (  2.0 );
  double signal  (  2.5 );
  parse_arguments(argc, argv, nRounds, alpha, scale, omega, signal);

  std::clog << "MAIN: settings are n=" << nRounds << " scale=" << scale << " omega=" << omega << " mu=" << signal << std::endl;
    
  // set up utility and wealth; angle 0.0 is not used
  WealthArray bidderWealth(omega, omega, nRounds, ScaledUniversalDist(scale));
  RiskVectorUtility utility(0.0, alpha);

  // open output file if want further output, otherwise direct to stdio
  if (false) {
    std::ostringstream ss;
    ss << "calc_risk_" << nRounds << "_" << signal << "_" << scale << "_" << omega ;
    std::ios_base::openmode mode = std::ios_base::trunc;
    std::ofstream output (ss.str(), mode);
  }
  // find expected risk, varying p_zero; write to output file
  for (double p=0.01; p < 1; p+=0.01)
  { std::pair<double,double> result (find_process_risk (nRounds, p, signal, utility, bidderWealth));
    std::cout << p << " " << result.first << " " << result.second << std::endl;
    // output << result.first << " " << result.second << std::endl;
  }
  
  return 0;
}


///////////////////



void
parse_arguments(int argc, char** argv,
		int &nRounds,  double &alpha,
		double &scale, double &omega,
		double &signal)
{
  static struct option long_options[] = {
    {"alpha",      required_argument, 0, 'A'},
    {"scale",      required_argument, 0, 's'},
    {"rounds",     required_argument, 0, 'n'},
    {"signal",     required_argument, 0, 'S'},
    {"omega",      required_argument, 0, 'W'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "A:s:n:N:p:S:W:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'A' : 
      {
	alpha = read_utils::lexical_cast<double>(optarg);
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

