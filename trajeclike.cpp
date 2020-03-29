#include <iostream>
#include <boost/program_options.hpp>
#include "trajec.h"

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
    ("maxlen,L", po::value<int>()->default_value(10), "max length of chop zone");

    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc,argv).options(opts).run();
    po::store (parsed, vm);
    po::notify(vm);

    if (vm.count("help")) {
      cout << opts << endl;
      return EXIT_SUCCESS;
    }

    const IndelParams params (vm.at("gamma").as<double>(), vm.at("mu").as<double>(), vm.at("r").as<double>());
    const ChopZoneConfig config (vm.at("maxevents").as<int>(), vm.at("maxlen").as<int>(), vm.at("verbose").as<int>());
    const double time = vm.at("time").as<double>();
    
    auto probs = chopZoneLikelihoods (params, time, config);
    double total = 0;
    for (const auto& pd: probs) {
      for (double p: pd)
	total += p;
      cout << to_string_join (pd) << endl;
    }
    if (config.verbose) {
      cerr << "Entry in row i, column j is probability of deleting i residues and inserting j residues before the next match" << endl;
      cerr << "Total: " << total << endl;
    }

    return EXIT_SUCCESS;
}
