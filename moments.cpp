#include <cmath>
#include <iostream>
#include "moments.h"
#include "rk4.h"
#include "util.h"

namespace TrajectoryLikelihood {

  Moments::Moments (const IndelParams& params, double tMax, double dt, int verbose) :
    params (params),
    mu (params.totalRightwardDeletionRatePerSite()),
    lambda (params.totalInsertionRatePerSite()),
    x (params.gamma * params.rIns),
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
    vector<vector<double> > pGap (maxLen + 1, vector<double> (maxLen + 1, 0));
    const auto& T = this->T.back();
    const double t = tMax;
    const double a = this->a(T,t), b = this->b(T,t), c = this->c(T,t);
    const double p = this->p(T,t), q = this->q(T,t), r = this->r(T,t);
    const double f = this->f(T,t), g = this->g(T,t), h = this->h(T,t);
    if (verbose > 2)
      cerr << "Pair HMM probability parameters: a=" << a << " b=" << b << " c=" << c
	   << " f=" << f << " g=" << g << " h=" << h
	   << " p=" << p << " q=" << q << " r=" << r
	   << endl;
    vector<vector<vector<double> > > fwd (2, vector<vector<double> > (maxLen + 1, vector<double> (2, 0)));  // fwd[i%2][j][state] where i=#ins, j=#del, state = 0(ins) or 1(del)
    for (int i = 0; i <= maxLen; ++i) {
      const int row = i % 2;
      const int prev = (i - 1) % 2;
      for (int j = 0; j <= maxLen; ++j) {
	fwd[row][j][0] = (i == 0
			  ? 0
			  : ((i == 1 && j == 0 ? b : 0)
			     + fwd[prev][j][1] * q
			     + fwd[prev][j][0] * g));

	fwd[row][j][1] = (j == 0
			  ? 0
			  : ((j == 1 && i == 0 ? c : 0)
			     + fwd[row][j-1][0] * h
			     + fwd[row][j-1][1] * r));
	
	pGap[i][j] = (i == 0 && j == 0
		      ? a
		      : (fwd[row][j][0] * f + fwd[row][j][1] * p));
      }
    }
    return pGap;
  }
}
