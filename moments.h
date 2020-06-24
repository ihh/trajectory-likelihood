#ifndef MOMENTS_INCLUDED
#define MOMENTS_INCLUDED

#include <vector>
#include "trajec.h"

namespace TrajectoryLikelihood {

  using namespace std;
  
  struct Moments {
    const IndelParams& params;
    const double mu, e;
    vector<vector<double> > MC;
    double dt, tMax;
    int verbose;
    Moments (const IndelParams&, double tMax, double dt, int verbose);
    double L (double t) const;
    double p (double M, double C, double t) const;
    double q (double M, double C, double t) const;
    double r (double M, double C, double t) const;
    double alpha (double M, double C, double t) const;
    double beta (double M, double C, double t) const;

    vector<vector<double> > chopZoneLikelihoods (int maxLen) const;
  };
}

#endif /* MOMENTS_INCLUDED */
