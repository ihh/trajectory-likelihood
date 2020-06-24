#include <cmath>
#include <iostream>
#include "moments.h"
#include "rk4.h"

namespace TrajectoryLikelihood {

  Moments::Moments (const IndelParams& params, double tMax, double dt, int verbose) :
    params (params),
    mu (params.totalRightwardDeletionRatePerSite()),
    e (params.r),
    MC (tMax / dt + 1, vector<double> (2)),
    dt (dt),
    tMax (tMax),
    verbose (verbose)
  {
    if (params.gamma != 1)
      throw runtime_error ("must have gamma=1 for moment-based approach");

    if (verbose > 1)
      cerr << "Moments method: rate=" << mu << " e=" << e << endl;

    auto eval_dmc_dt = [&] (double t, const vector<double>& mc, vector<double>& dmc_dt) {
      const double M = mc[0];
      const double C = mc[1];
      const double q_t = q (M, C, t);
      const double r_t = r (M, C, t);
      const double alpha_t = alpha (M, C, t);
      const double beta_t = beta (M, C, t);
      dmc_dt[0] = mu * (-2*M + (1-M)*(1-e)*(1-q_t)*(1-r_t)/(2*(1-e*q_t)));
      dmc_dt[1] = mu * alpha_t / beta_t;
    };

    MC[0][0] = 1.;
    MC[0][1] = 0.;
    RungeKutta4 mc_solver (2);
    for (int i = 1; i < MC.size(); ++i) {
      const double t = dt * (double) i;
      mc_solver.step (eval_dmc_dt, t, dt, MC[i-1], MC[i]);
      if (verbose > 1)
	cerr << "t=" << t << " M=" << MC[i][0] << " L=" << L(t) << " C=" << MC[i][1] << endl;
    }
  }

  double Moments::L (double t) const {
    return exp (mu*t/(1-e)) - 1;
  }
  
  double Moments::p (double M, double C, double t) const {
    return M;
  }
  
  double Moments::q (double M, double C, double t) const {
    const double L = this->L(t);
    return (4*L*L*(2*L+M-1) + C*(M-1)*(4*L+M-1)) / (8*L*L*L - 4*C*L*(M-1));
  }
  
  double Moments::r (double M, double C, double t) const {
    const double L = this->L(t);
    return (1-M)*C/(4*L*L - (1-M)*C);
  }
  
  double Moments::alpha (double M, double C, double t) const {
    const double L = this->L(t);
    return C * (1 - e) * (C * (M - 1) + L * (3 * (M - 1) + 2 * L * (1 + M))) + 2 * L * L * (2 * L * (3 + 2 * L) - e * (L * M + L + M - 1));
  }
  
  double Moments::beta (double M, double C, double t) const {
    const double L = this->L(t);
    return (1 - e) * (C * (1 - e) * (M - 1) + 2 * L * (2 * L - e * (L * M + L + M - 1)));
  }

  vector<vector<double> > Moments::chopZoneLikelihoods (int maxLen) const {
    vector<vector<double> > pGap (maxLen + 1, vector<double> (maxLen + 1, 0));
    const auto& mc = MC.back();
    const double M = mc[0], C = mc[1];
    const double t = tMax;
    const double p = this->p(M,C,t), q = this->q(M,C,t), r = this->r(M,C,t);
    if (verbose)
      cerr << "p=" << p << " q=" << q << " r=" << r << endl;
    vector<vector<vector<double> > > fwd (2, vector<vector<double> > (maxLen + 1, vector<double> (2, 0)));  // fwd[i%2][j][state] where i=#ins, j=#del, state = 0(ins) or 1(del)
    for (int i = 0; i <= maxLen; ++i) {
      const int row = i % 2;
      const int prev = (i - 1) % 2;
      for (int j = 0; j <= maxLen; ++j) {
	fwd[row][j][0] = (i == 0
			  ? 0
			  : ((i == 1 && j == 0 ? ((1 - p) / 2) : 0)
			     + fwd[prev][j][1] * (1-q) * r
			     + fwd[prev][j][0] * q));

	fwd[row][j][1] = (j == 0
			  ? 0
			  : ((j == 1 && i == 0 ? ((1 - p) / 2) : 0)
			     + fwd[row][j-1][0] * (1-q) * r
			     + fwd[row][j-1][1] * q));
	
	pGap[i][j] = (i == 0 && j == 0
		      ? p
		      : ((fwd[row][j][0] + fwd[row][j][1]) * (1-q) * (1-r)));
      }
    }
    return pGap;
  }
}
