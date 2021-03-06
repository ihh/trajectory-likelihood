#include <chrono>
#include <iostream>
#include <boost/program_options.hpp>
#include "trajec.h"
#include "simulate.h"
#include "moments.h"
#include "cim.h"
#include "tkf.h"

using namespace TrajectoryLikelihood;
using namespace std;
namespace po = boost::program_options;

int main (int argc, char** argv) {
  po::options_description opts ("Options");
  opts.add_options()
    ("help,h", "display this help message")
    ("verbose,v", po::value<int>()->default_value(0), "logging verbosity")
    ("lambda,L", po::value<double>(), "insertion rate (default is same as mu)")
    ("mu,M", po::value<double>()->default_value(.049), "deletion rate")
    ("x,X", po::value<double>(), "insertion extension probability (default is same as y)")
    ("y,Y", po::value<double>()->default_value(.543), "deletion extension probability")
    ("time,t", po::value<double>()->default_value(1), "time parameter")
    ("maxevents,E", po::value<int>()->default_value(3), "max # of indel events in trajectory")
    ("maxlen,l", po::value<int>()->default_value(10), "max length of chop zone")
    ("benchmark,b", "perform time benchmark instead of likelihood calculation")
    ("simulate,s", "perform stochastic simulation instead of likelihood calculation")
    ("counts,c", "report simulation counts instead of probabilities")
    ("initlen,i", po::value<int>()->default_value(1000), "initial sequence length for simulation")
    ("trials,n", po::value<int>()->default_value(100000), "number of simulation trials")
    ("moments,m", "use method of moments (Holmes 2020) for likelihood calculations")
    ("cim,C", "use De Maio 2020 Cumulative Indel Model for likelihood calculations")
    ("tkf91", "use TKF91 Pair HMM for likelihood calculations")
    ("tkf92", "use TKF92 Pair HMM for likelihood calculations")
    ("rs07", "use RS07 Pair HMM for likelihood calculations")
    ("prank", "use PRANK Pair HMM for likelihood calculations")
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

    // logging
    const int verbose = vm.at("verbose").as<int>();

    // get General Geometric Indel Model params (De Maio 2020; Holmes 2020)
    const double mu_ggi = vm.at("mu").as<double>();
    const double lambda_ggi = vm.count("lambda") ? vm.at("lambda").as<double>() : mu_ggi;
    const double rDel = vm.at("y").as<double>();
    const double rIns = vm.count("x") ? vm.at("x").as<double>() : rDel;

    // convert params to Long Indel Model scaling (Miklos, Lunter & Holmes 2004)
    const double mu_li = mu_ggi / (1. - rDel);
    const double gamma_li = (lambda_ggi * (1. - rIns)) / (mu_ggi * (1. - rDel));
    const IndelParams params (gamma_li, mu_li, rDel, rIns);

    if (verbose) {
      // log params in GGI format
      cerr << "mu_GGI=" << params.totalRightwardDeletionRatePerSite() << " lambda_GGI=" << params.totalInsertionRatePerSite() << " x=" << params.rIns << " y=" << params.rDel << endl;
    }

    // get remaining params
    const double t = vm.at("time").as<double>();

    const bool reportCounts = vm.count("simulate") && vm.count("counts");
    const int maxLen = vm.at("maxlen").as<int>();
    const double dt = vm.at("dt").as<double>();
    vector<vector<double> > probs;

    if (vm.count("benchmark")) {
      const int E = vm.at("maxevents").as<int>();
      cout << "chopZoneLength meanTime/uS error" << endl;
      const int trials = vm.at("trials").as<int>();
      for (int L = 1; L <= maxLen; ++L) {
	double tSum = 0, tSum2 = 0;
	for (int n = 0; n < trials; ++n) {
	  const auto tStart = std::chrono::system_clock::now();
	  const ChopZoneConfig config (E, L, verbose);
	  (void) chopZoneLikelihoods (params, t, config);
	  const auto tEnd = std::chrono::system_clock::now();
	  const double t = std::chrono::duration_cast< std::chrono::microseconds >(tEnd - tStart).count();
	  tSum += t;
	  tSum2 += t*t;
	}
	const double tMean = tSum / (double) trials, t2Mean = tSum2 / (double) trials;
	cout << L << " " << tMean << " " << sqrt ((t2Mean - tMean*tMean) / (double) trials) << endl;
      }
    } else {
    
      if (vm.count("simulate")) {
	const SimulationConfig config (vm.at("initlen").as<int>(), maxLen, vm.at("trials").as<int>(), verbose);
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
	const Moments moments (params, t, dt, verbose);
	probs = moments.chopZoneLikelihoods (maxLen);
      } else if (vm.count("cim")) {
	const CumulativeIndelModel cim (params, t, dt, verbose);
	probs = cim.chopZoneLikelihoods (maxLen);
      } else if (vm.count("tkf91")) {
	const TKF91 tkf91 (params, t, verbose);
	probs = tkf91.chopZoneLikelihoods (maxLen);
      } else if (vm.count("tkf92")) {
	const TKF92 tkf92 (params, t, verbose);
	probs = tkf92.chopZoneLikelihoods (maxLen);
      } else if (vm.count("prank")) {
	const PRANK prank (params, t, verbose);
	probs = prank.chopZoneLikelihoods (maxLen);
      } else if (vm.count("rs07")) {
	const RS07 rs07 (params, t, verbose);
	probs = rs07.chopZoneLikelihoods (maxLen);
      } else {
	const ChopZoneConfig config (vm.at("maxevents").as<int>(), maxLen, verbose);
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
