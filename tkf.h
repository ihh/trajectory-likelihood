#ifndef TKF_INCLUDED
#define TKF_INCLUDED

#include <vector>
#include "trajec.h"

namespace TrajectoryLikelihood {

  using namespace std;

  struct TKF {
    const IndelParams& params;
    int verbose;
    const double mu, lambda, alpha, beta, gamma;
    TKF (const IndelParams&, double t, int verbose);
  };

  struct TKF91 : TKF {
    TKF91 (const IndelParams&, double t, int verbose);
    vector<vector<double> > chopZoneLikelihoods (int maxLen) const;
  };

  struct TKF92 : TKF {
    double r;
    TKF92 (const IndelParams&, double t, int verbose);
    vector<vector<double> > chopZoneLikelihoods (int maxLen) const;
  };

  struct PRANK {
    int verbose;
    double epsilon, gamma, delta;
    PRANK (const IndelParams&, double t, int verbose);
    vector<vector<double> > chopZoneLikelihoods (int maxLen) const;
  };

  struct RS07 {
    int verbose;
    double epsilon, delta;
    RS07 (const IndelParams&, double t, int verbose);
    vector<vector<double> > chopZoneLikelihoods (int maxLen) const;
  };
}

#endif /* TKF_INCLUDED */
