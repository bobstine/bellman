#include "bellman.h"
#include "random.h"
#include "read_utils.h"  // coerce

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <ios>

#include <algorithm> // generate


inline std::ostream&
operator<<(std::ostream &output, const std::pair<double,double> pair)
{ output << " < " << pair.first << " , " << pair.second << " > "; return output; }




void
parse_arguments(int argc, char** argv,
		int &nRounds, int &nReplications, double &angle, double &alpha,
		double &scale, double &omega,double &pZero, double &signal);


// Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     Main     
int  main(int argc, char** argv)
{
  // set default control values
  int    nRounds (100   );
  int    nReps   ( 25   );
  double omega   (  0.5 );
  double alpha   (  0.05);
  double scale   (  2.0 );
  double angle   (210   );
  double signal  (  2.5 );
  double pZero   (  0.50);
  parse_arguments(argc, argv, nRounds, nReps, angle, alpha, scale, omega, pZero, signal);

  std::clog << "MAIN: settings are n=" << nRounds << "  N=" << nReps << " theta=" << angle << " scale=" << scale
	    << "  pZero=" << pZero << " mu=" << signal << std::endl;
    
    
  // set up utility and wealth
  WealthArray bidderWealth(omega, omega, nRounds, ScaledUniversalDist(scale));
  RiskVectorUtility utility(angle, alpha);

  // generate means
  RandomGenerator randu(23984);
  std::vector<double> means (nRounds);

  // open output file
  std::ostringstream ss;
  ss << "sim_risk_" << angle << "_" << pZero << "_" << signal;
  std::ios_base::openmode mode = std::ios_base::trunc;
  std::ofstream output (ss.str(), mode);

  // run simulation
  for (int rep = 0; rep < nReps; ++rep)
  { std::generate(means.begin(), means.end(), [&randu, pZero, signal]()->double { if(randu.uniform() < pZero) return 0.0; else return signal; });
    std::pair<double,double> result (find_process_risk (utility, bidderWealth, means));
    if (0 == (rep%10)) std::clog << rep << " ";
    output << result.first << " " << result.second << std::endl;
  }
  return 0;
}


///////////////////



void
parse_arguments(int argc, char** argv,
		int &nRounds, int &nReps,
		double &angle, double &alpha,
		double &scale, double &omega,
		double &pZero, double &signal)
{
  static struct option long_options[] = {
    {"angle",      required_argument, 0, 'a'},
    {"alpha",      required_argument, 0, 'A'},
    {"scale",      required_argument, 0, 's'},
    {"length",     required_argument, 0, 'n'},
    {"reps",       required_argument, 0, 'N'},
    {"pzero",      required_argument, 0, 'p'},
    {"signal",     required_argument, 0, 'S'},
    {"omega",      required_argument, 0, 'W'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "a:A:s:n:N:p:S:W:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'a' : 
      {
	angle = read_utils::lexical_cast<double>(optarg);
	break;
      }
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
    case 'N' :
      {
	nReps = read_utils::lexical_cast<int>(optarg);
	break;
      }
    case 'p' : 
      {
	pZero = read_utils::lexical_cast<double>(optarg);
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

