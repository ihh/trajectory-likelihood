#include <cmath>
#include <iostream>
#include "moments.h"
#include "rk4.h"
#include "util.h"
#include "mid.h"

namespace TrajectoryLikelihood {

  Moments::Moments (const IndelParams& params, double tMax, double dt, int verbose) :
    params (params),
    mu (params.totalRightwardDeletionRatePerSite()),
    lambda (params.totalInsertionRatePerSite()),
    x (params.rIns),
    y (params.rDel),
    T (tMax / dt + 1, vector<double> (4, 0.)),
    dt (dt),
    tMax (tMax),
    verbose (verbose)
  {
    auto eval_dmc_dt = [&] (double t, const vector<double>& T, vector<double>& dT_dt) {
      const double a = this->a(T,t), b = this->b(T,t), c = this->c(T,t);
      const double p = this->p(T,t), q = this->q(T,t), r = this->r(T,t);
      const double f = this->f(T,t), g = this->g(T,t), h = this->h(T,t);
      dT_dt[0] = (-1 + b + c) * (lambda + mu) - (b * f * mu * (-1 + y))/(1 + (-1 + f + h) * y);
      dT_dt[1] = lambda - b * lambda - (b * (f + h) * mu)/(1 + (-1 + f + h) * y);
      dT_dt[2] = -(-1 + b + c) * lambda - (f * (f + h) * mu * (c * q + b * (p + q)))/((h * p + f * (p + q)) * (1 + (-1 + f + h) * y));
      dT_dt[3] = ((f + h) * mu * (-c * (-1 + f + h) * q + b * (p + q - h * q)))/((h * p + f * (p + q)) * (1 + (-1 + f + h) * y));
      if (verbose > 3)
	cerr << "t=" << t << " T=" << vec_to_string(T)
	     << " SI=" << SI(t) << " SD=" << SD(t)
	     << " a=" << a << " b=" << b << " c=" << c
	     << " f=" << f << " g=" << g << " h=" << h
	     << " p=" << p << " q=" << q << " r=" << r
	     << " dT/dt=" << vec_to_string(dT_dt)
	     << endl;
    };

    T[0][0] = 1.;
    RungeKutta4 T_solver (4);
    for (int i = 1; i < T.size(); ++i) {
      const double t = dt * (double) (i - 1);
      T_solver.step (eval_dmc_dt, t, dt, T[i-1], T[i]);
      if (verbose > 2)
	cerr << "t=" << t
	     << " T=" << vec_to_string(T[i])
	     << endl;
    }
  }

  vector<vector<double> > Moments::chopZoneLikelihoods (int maxLen) const {
    const auto& T = this->T.back();
    const double t = tMax;
    const double a = this->a(T,t), b = this->b(T,t), c = this->c(T,t);
    const double f = this->f(T,t), g = this->g(T,t), h = this->h(T,t);
    const double p = this->p(T,t), q = this->q(T,t), r = this->r(T,t);
    const MID_HMM hmm (a, b, c, f, g, h, p, q, r, verbose);
    return hmm.chopZoneLikelihoods (maxLen);
  }
}
