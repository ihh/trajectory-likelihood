#include <cmath>
#include <iostream>
#include "cim.h"
#include "rk4.h"
#include "mid.h"

namespace TrajectoryLikelihood {

  CumulativeIndelModel::CumulativeIndelModel (const IndelParams& params, double tMax, double dt, int verbose) :
    params (params),
    ri (params.totalInsertionRatePerSite()),
    rd (params.totalRightwardDeletionRatePerSite()),
    gi (params.rIns),
    gd (params.rDel),
    A (tMax / dt + 1, vector<double> (3, 0.)),
    dt (dt),
    tMax (tMax),
    verbose (verbose)
  {
    if (verbose > 1)
      cerr << "Cumulative indel method: ri=" << ri << " gi=" << gi << endl;

    auto eval_dA_dt = [&] (double t, const vector<double>& A, vector<double>& dA_dt) {
      const double Ai = A[0];
      const double Ad = A[1];
      const double Aid = A[2];
      const double Pm = Pmt(t);
      const double Pi = Pit(A,t);
      const double Pid = Pidt(A,t);
      const double _git = git(A,t);
      const double Li = Lit(t);
      const double denom = (1 - gd * (1 - Pi)) * (1 - gd * _git);
      dA_dt[0] = (Pm - Ai) * ri - Ai * rd / (1 - gd) + (Pm - Ai) * rd * Pi * (1 - gd) / denom - Ai * rd * (1 - _git) * (1 - gd) / denom;
      dA_dt[1] = (Pm - Ad) * rd - Ad * rd / (1 - gd);
      if (Pi != 0)
	dA_dt[1] += ((Pi - Pid) / Pi) * (Li * rd * gd * (1 - _git)) / (1 - gd * _git);
      dA_dt[2] = (Ad - Aid) * ri - Aid * rd / (1 - gd) + (Ai - Aid) * rd / (1 - _git * gd) + (Pm - Ai) * rd * Pi * (1 - gd) / denom - Aid * rd * (1 - _git) * (1 - gd) / denom + (Ai - Aid) * rd * (Pi * (1 - gd) / denom) * gd * (1 - _git) / (1 - gd * _git);
    };

    RungeKutta4 A_solver (3);
    for (int i = 1; i < A.size(); ++i) {
      const double t = dt * (double) i;
      A_solver.step (eval_dA_dt, t, dt, A[i-1], A[i]);
      if (verbose > 2)
	cerr << "t=" << t << " Ai=" << A[i][0] << " Ad=" << A[i][1] << " Aid=" << A[i][2] << endl;
    }
  }

  double CumulativeIndelModel::Lit (double t) const {
    return Pmt(t) * (exp (ri*t/(1-gi)) - 1);
  }
  
  double CumulativeIndelModel::Pmt (double t) const {
    return exp (-t*rd/(1-gd));
  }
  
  double CumulativeIndelModel::Pit (const vector<double>& A, double t) const {
    return A[0] / Pmt(t);
  }
  
  double CumulativeIndelModel::Pdt (const vector<double>& A, double t) const {
    return A[1] / Pmt(t);
  }
  
  double CumulativeIndelModel::Pidt (const vector<double>& A, double t) const {
    return A[2] / Pmt(t);
  }

  double CumulativeIndelModel::git (const vector<double>& A, double t) const {
    return 1 - A[0] / Lit(t);
  }
  
  double CumulativeIndelModel::gdt (const vector<double>& A, double t) const {
    return 1 - A[1] / (1 - Pmt(t));
  }

  vector<vector<double> > CumulativeIndelModel::chopZoneLikelihoods (int maxLen) const {
    const auto& a = A.back();
    const double t = tMax;
    const double m2m = 1 - Pit(a,t) - Pdt(a,t) + Pidt(a,t),
      m2d = Pdt(a,t),
      m2i = Pit(a,t) - Pidt(a,t),
      d2m = (1 - gdt(a,t)) * (Pdt(a,t) - Pidt(a,t)) / Pdt(a,t),
      d2d = gdt(a,t),
      d2i = (1 - gdt(a,t)) * Pidt(a,t) / Pdt(a,t),
      i2m = 1 - git(a,t),
      i2i = git(a,t);
    const MID_HMM hmm (m2m, m2i, m2d, i2m, i2i, 0., d2m, d2i, d2d, verbose);
    return hmm.chopZoneLikelihoods (maxLen);
  }
}
