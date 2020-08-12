#ifndef MID_INCLUDED
#define MID_INCLUDED

#include <vector>

namespace TrajectoryLikelihood {

  using namespace std;

  struct MID_HMM {
    const double a, b, c, f, g, h, p, q, r;
    MID_HMM (double m2m, double m2i, double m2d, double i2m, double i2i, double i2d, double d2m, double d2i, double d2d, int verbose = 0);
    vector<vector<double> > chopZoneLikelihoods (int maxLen) const;
  };
}

#endif /* MID_INCLUDED */
