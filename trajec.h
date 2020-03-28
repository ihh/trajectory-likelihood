#include <machineboss/machine.h>
#include <machineboss/forward.h>
#include <machineboss/util.h>

using namespace std;
using namespace MachineBoss;

// Algorithm 1 of (Miklos, Lunter & Holmes, 2004)
// exitRates = chi
// transitionRates = r
double trajectoryLikelihood (const vector<double>& exitRates,
			     const vector<double>& transitionRates,
			     double time);

// Parameters
struct IndelParams {
  double gamma, mu, r;
  IndelParams() : gamma(1), mu(0), r(0) { }
  IndelParams (double g, double m, double r) : gamma(g), mu(m), r(r) { }
};
  
// Config for chop zone probability calculations
struct ChopZoneConfig {
  int maxEvents;
  int maxLen;
  bool verbose;
  ChopZoneConfig() : maxEvents(3), maxLen(100), verbose(false) { }
  ChopZoneConfig (int e, int l, bool v) : maxEvents(e), maxLen(l), verbose(v) { }
};

// Calculate chop zone probabilities (internal zones i.e. not at the ends of the sequence)
double chopZoneLikelihood (int nDeleted, int nInserted, const IndelParams& params, double time, const ChopZoneConfig& config);
vector<vector<double> > chopZoneLikelihoods (const IndelParams& params, double time, const ChopZoneConfig& config);
