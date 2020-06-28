#include <chrono>
#include <iostream>
#include <boost/program_options.hpp>
#include "trajec.h"
#include "simulate.h"
#include "moments.h"
#include "cim.h"

using namespace TrajectoryLikelihood;
using namespace std;
namespace po = boost::program_options;

int main (int argc, char** argv) {
  po::options_description opts ("Options");
  opts.add_options()
    ("help,h", "display this help message")
    ("verbose,v", po::value<int>()->default_value(0), "logging verbosity")
    ("gamma", po::value<double>()->default_value(.99), "gamma parameter from MLH 2004 (parameter of equilibrium geometric distribution over sequence lengths)")
    ("mu", po::value<double>()->default_value(.049), "mu parameter from MLH 2004 (rightward deletion rate)")
    ("r,r", po::value<double>()->default_value(.543), "r parameter from MLH 2004 (parameter of geometric distribution over deletion lengths)")
    ("rins,I", po::value<double>(), "(1/gamma) * parameter of geometric distribution over insertion lengths (default is same as deletion rate)")
    ("time,t", po::value<double>()->default_value(1), "time parameter")
    ("maxevents,E", po::value<int>()->default_value(3), "max # of indel events in trajectory")
    ("maxlen,L", po::value<int>()->default_value(10), "max length of chop zone")
    ("benchmark,b", "perform time benchmark instead of likelihood calculation")
    ("simulate,s", "perform stochastic simulation instead of likelihood calculation")
    ("counts,c", "report simulation counts instead of probabilities")
    ("initlen,i", po::value<int>()->default_value(1000), "initial sequence length for simulation")
    ("trials,n", po::value<int>()->default_value(100000), "number of simulation trials")
    ("moments,m", "use method of moments for likelihood calculations")
    ("cim,C", "use de Maio's Cumulative Indel Model for likelihood calculations")
    ("dt,D", po::value<double>()->default_value(.01), "time step for numerical integration")
    ("seed,d", po::value<unsigned long long>()->default_value(mt19937::default_seed), "seed for random number generator")
    ("seedtime,T", "use current time as a seed");

    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc,argv).options(opts).run();
    po::store (parsed, vm);
    po::notify(vm);

    if (vm.count("help")) {
      cout << opts << endl;
      return EXIT_SUCCESS;
    }

    const double gamma = vm.at("gamma").as<double>();
    const double mu = vm.at("mu").as<double>();
    const double rDel = vm.at("r").as<double>();
    const IndelParams params (gamma, mu, rDel, vm.count("rins") ? vm.at("rins").as<double>() : rDel);
    const double t = vm.at("time").as<double>();

    const int verbose = vm.at("verbose").as<int>();

    const bool reportCounts = vm.count("simulate") && vm.count("counts");
    vector<vector<double> > probs;

    if (vm.count("benchmark")) {
      const int E = vm.at("maxevents").as<int>(), maxL = vm.at("maxlen").as<int>();
      if (verbose)
	cerr << "chopZoneLength time/uS" << endl;
      for (int L = 1; L <= maxL; ++L) {
	const auto tStart = std::chrono::system_clock::now();
	const ChopZoneConfig config (E, L, verbose);
	(void) chopZoneLikelihoods (params, t, config);
	const auto tEnd = std::chrono::system_clock::now();
	cout << L << " " << std::chrono::duration_cast< std::chrono::microseconds >(tEnd - tStart).count() << endl;
      }
    } else {
    
      if (vm.count("simulate")) {
	const SimulationConfig config (vm.at("initlen").as<int>(), vm.at("maxlen").as<int>(), vm.at("trials").as<int>(), verbose);
	const bool useClock = vm.count("seedtime");
	unsigned long long seed;
	if (useClock)
	  seed = std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch()).count();
	else
	  seed = vm.at("seed").as<unsigned long long>();
	cerr << "Random number seed: " << seed << endl;
	mt19937 rnd (seed);
	probs = chopZoneSimulatedProbabilities (params, t, config, rnd, reportCounts);
      } else if (vm.count("moments")) {
	const Moments moments (params, t, vm.at("dt").as<double>(), verbose);
	probs = moments.chopZoneLikelihoods (vm.at("maxlen").as<int>());
      } else if (vm.count("cim")) {
	const CumulativeIndelModel cim (params, t, vm.at("dt").as<double>(), verbose);
	probs = cim.chopZoneLikelihoods (vm.at("maxlen").as<int>());
      } else {
	const ChopZoneConfig config (vm.at("maxevents").as<int>(), vm.at("maxlen").as<int>(), verbose);
	probs = chopZoneLikelihoods (params, t, config);
      }

      double total = 0;
      for (const auto& pd: probs) {
	for (double p: pd)
	  total += p;
	cout << to_string_join (pd) << endl;
      }
      if (verbose) {
	cerr << "Entry in row i, column j is " << (reportCounts ? "count" : "probability") << " of inserting i residues and deleting j residues before the next match" << endl;
	if (reportCounts)
	  cerr << "Final count represents overflow (chop zones that were larger than size limit)" << endl;
	cerr << "Total: " << total << endl;
	if (!vm.count("counts")) {
	  double ei = 0, ed = 0, ei2 = 0, ed2 = 0, eid = 0, pi0 = 0, pd0 = 0;
	  for (int i = 0; i < probs.size(); ++i)
	    for (int d = 0; d < probs[i].size(); ++d) {
	      const double p = probs[i][d];
	      ei += i * p;
	      ed += d * p;
	      ei2 += i * i * p;
	      ed2 += d * d * p;
	      eid += i * d * p;
	      if (i == 0)
		pi0 += p;
	      if (d == 0)
		pd0 += p;
	    }
	  cerr << "E[#ins]=" << ei << " E[#del]=" << ed
	       << " E[#ins^2]=" << ei2 << " E[#del^2]=" << ed2
	       << " V[#ins]=" << (ei2 - ei*ei) << " V[#del]=" << (ed2 - ed*ed)
	       << " E[#ins*#del]=" << eid
	       << " Cov(#ins,#del)=" << (eid-ei*ed)
	       << endl;
	}
      }
    }
    
    return EXIT_SUCCESS;
}
