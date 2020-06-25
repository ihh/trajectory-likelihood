#include <iostream>
#include <boost/program_options.hpp>
#include "trajec.h"
#include "simulate.h"
#include "moments.h"

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
    ("r,r", po::value<double>()->default_value(.543), "r parameter from MLH 2004 (parameter of geometric distribution over deletion lengths")
    ("time,t", po::value<double>()->default_value(1), "time parameter")
    ("maxevents,E", po::value<int>()->default_value(3), "max # of indel events in trajectory")
    ("maxlen,L", po::value<int>()->default_value(10), "max length of chop zone")
    ("simulate,s", "perform stochastic simulation instead of likelihood calculation")
    ("counts,c", "report simulation counts instead of probabilities")
    ("initlen,i", po::value<int>()->default_value(1000), "initial sequence length for simulation")
    ("trials,n", po::value<int>()->default_value(100000), "number of simulation trials")
    ("moments,m", "use method of moments for likelihood calculations")
    ("dt,D", po::value<double>()->default_value(.01), "time step for numerical integration")
    ("seed,d", po::value<int>()->default_value(mt19937::default_seed), "seed for random number generator");

    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc,argv).options(opts).run();
    po::store (parsed, vm);
    po::notify(vm);

    if (vm.count("help")) {
      cout << opts << endl;
      return EXIT_SUCCESS;
    }

    const IndelParams params (vm.at("gamma").as<double>(), vm.at("mu").as<double>(), vm.at("r").as<double>());
    const double time = vm.at("time").as<double>();

    const int verbose = vm.at("verbose").as<int>();

    const bool reportCounts = vm.count("simulate") && vm.count("counts");
    vector<vector<double> > probs;

    if (vm.count("simulate")) {
      const SimulationConfig config (vm.at("initlen").as<int>(), vm.at("maxlen").as<int>(), vm.at("trials").as<int>(), verbose);
      mt19937 rnd (vm.at("seed").as<int>());
      probs = chopZoneSimulatedProbabilities (params, time, config, rnd, reportCounts);
    } else if (vm.count("moments")) {
      const Moments moments (params, time, vm.at("dt").as<double>(), verbose);
      probs = moments.chopZoneLikelihoods (vm.at("maxlen").as<int>());
    } else {
      const ChopZoneConfig config (vm.at("maxevents").as<int>(), vm.at("maxlen").as<int>(), verbose);
      probs = chopZoneLikelihoods (params, time, config);
    }

    double total = 0;
    for (const auto& pd: probs) {
      for (double p: pd)
	total += p;
      cout << to_string_join (pd) << endl;
    }
    if (verbose) {
      cerr << "Entry in row i, column j is " << (reportCounts ? "count" : "probability") << " of deleting i residues and inserting j residues before the next match" << endl;
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

    return EXIT_SUCCESS;
}
