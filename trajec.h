#include <vector>

namespace TrajectoryLikelihood {

  using namespace std;

  // Algorithm 1 of (Miklos, Lunter & Holmes, 2004)
  // Calculates the likelihood of a state trajectory in a continuous-time Markov chain, integrating out the event times.
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
    int verbose;
    ChopZoneConfig() : maxEvents(3), maxLen(100), verbose(0) { }
    ChopZoneConfig (int e, int l, int v) : maxEvents(e), maxLen(l), verbose(v) { }
  };

  // Calculate chop zone probabilities, i.e. gap probabilities in the long indel model
  // Currently only implemented for internal zones i.e. not the zones at the ends of the sequence.
  double chopZoneLikelihood (int nDeleted, int nInserted, const IndelParams& params, double time, const ChopZoneConfig& config);
  vector<vector<double> > chopZoneLikelihoods (const IndelParams& params, double time, const ChopZoneConfig& config);   // result is indexed [nDeleted][nInserted]

  // Templates
  // to_string_join
  template<class Container>
  std::string to_string_join (const Container& c, const char* sep = " ") {
    std::ostringstream j;
    int n = 0;
    for (const auto& s : c) {
      if (n++ > 0)
	j << sep;
      j << s;
    }
    return j.str();
  }

}
