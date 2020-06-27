#ifndef MOMENTS_INCLUDED
#define MOMENTS_INCLUDED

#include <vector>
#include "trajec.h"

namespace TrajectoryLikelihood {

  using namespace std;
  
  struct Moments {
    const IndelParams& params;
    const double mu, lambda, x, y;
    vector<vector<double> > T;  // TMM, TMI, TIM, TDI
    double dt, tMax;
    int verbose;
    Moments (const IndelParams&, double tMax, double dt, int verbose);

    inline double SI (double t) const { return exp (lambda*t/(1-x)) - 1; }
    inline double SD (double t) const { return exp (mu*t/(1-y)) - 1; }

    inline double TMM (const vector<double>& T, double t) const { return T[0]; }
    inline double TMI (const vector<double>& T, double t) const { return T[1]; }
    inline double TMD (const vector<double>& T, double t) const { return 1 - T[0] - T[1]; }
    inline double TIM (const vector<double>& T, double t) const { return T[2]; }
    inline double TII (const vector<double>& T, double t) const { return SI(t) - TMI(T,t) - TDI(T,t); }
    inline double TID (const vector<double>& T, double t) const { return TMI(T,t) + TDI(T,t) - TIM(T,t); }
    inline double TDM (const vector<double>& T, double t) const { return 1 - TMM(T,t) - TIM(T,t); }
    inline double TDI (const vector<double>& T, double t) const { return T[3]; }
    inline double TDD (const vector<double>& T, double t) const { return SD(t) + TMM(T,t) + TIM(T,t) - TDI(T,t) - 1; }

    inline double a (const vector<double>& T, double t) const { return TMM(T,t); }
    inline double b (const vector<double>& T, double t) const { return TMI(T,t); }
    inline double c (const vector<double>& T, double t) const { return TMD(T,t); }
    inline double f (const vector<double>& T, double t) const { return TIM(T,t) / SI(t); }
    inline double g (const vector<double>& T, double t) const { return TII(T,t) / SI(t); }
    inline double h (const vector<double>& T, double t) const { return TID(T,t) / SI(t); }
    inline double p (const vector<double>& T, double t) const { return TDM(T,t) / SD(t); }
    inline double q (const vector<double>& T, double t) const { return TDI(T,t) / SD(t); }
    inline double r (const vector<double>& T, double t) const { return TDD(T,t) / SD(t); }
    
    vector<vector<double> > chopZoneLikelihoods (int maxLen) const;
  };
}

#endif /* MOMENTS_INCLUDED */
