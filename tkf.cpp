#include <cmath>
#include <iostream>
#include "tkf.h"
#include "mid.h"

namespace TrajectoryLikelihood {

  TKF::TKF (const IndelParams& params, double t, int verbose) :
    params (params), verbose (verbose),
    mu (params.mu),
    lambda (params.gamma * params.mu),
    alpha (exp(-mu*t)),
    beta ((lambda*(exp(-lambda*t)-exp(-mu*t)))/(mu*exp(-lambda*t)-lambda*exp(-mu*t))),
    gamma (1 - mu*beta/(lambda*(1-alpha)))
  {
    if (verbose > 4)
      cerr << "TKF rates: lambda=" << lambda << " mu=" << mu << endl
	   << "TKF probabilities: alpha=" << alpha << " beta=" << beta << " gamma=" << gamma << endl;
  }

  TKF91::TKF91 (const IndelParams& params, double t, int verbose) :
    TKF (params, t, verbose)
  { }

  vector<vector<double> > TKF91::chopZoneLikelihoods (int maxLen) const {
    const MID_HMM hmm ((1-beta)*alpha, beta, (1-beta)*(1-alpha),
		       (1-beta)*alpha, beta, (1-beta)*(1-alpha),
		       (1-gamma)*alpha, gamma, (1-gamma)*(1-alpha),
		       verbose);
    return hmm.chopZoneLikelihoods (maxLen);
  }

  TKF92::TKF92 (const IndelParams& params, double t, int verbose) :
    TKF (params, t, verbose),
    r (params.rDel)
  { }

  vector<vector<double> > TKF92::chopZoneLikelihoods (int maxLen) const {
    const MID_HMM hmm (r + (1-r)*(1-beta)*alpha, (1-r)*beta, (1-r)*(1-beta)*(1-alpha),
		       (1-r)*(1-beta)*alpha, r + (1-r)*beta, (1-r)*(1-beta)*(1-alpha),
		       (1-r)*(1-gamma)*alpha, (1-r)*gamma, r + (1-r)*(1-gamma)*(1-alpha),
		       verbose);
    return hmm.chopZoneLikelihoods (maxLen);
  }

  PRANK::PRANK (const IndelParams& params, double t, int verbose) :
    verbose (verbose),
    epsilon (params.rDel),
    gamma (epsilon),
    delta (1 - exp(-params.totalInsertionRatePerSite()*t/(1-gamma)))
  {
    if (verbose > 4)
      cerr << "PRANK probabilities: epsilon=" << epsilon << " gamma=" << gamma << " delta=" << delta << endl;
  }
  
  vector<vector<double> > PRANK::chopZoneLikelihoods (int maxLen) const {
    const MID_HMM hmm (gamma + (1-gamma)*(1-2*delta), (1-gamma)*delta, (1-gamma)*delta,
		       (1-epsilon)*(1-2*delta), epsilon + (1-epsilon)*delta, (1-epsilon)*delta,
		       (1-epsilon)*(1-2*delta), epsilon + (1-epsilon)*delta, (1-epsilon)*delta,
		       verbose);
    return hmm.chopZoneLikelihoods (maxLen);
    
  }

  RS07::RS07 (const IndelParams& params, double t, int verbose) :
    verbose (verbose),
    epsilon (params.rDel),
    delta (1 / (1 + 1 / (1 - exp(-params.totalInsertionRatePerSite()*t/(1-epsilon)))))
  {
    if (verbose > 4)
      cerr << "PRANK probabilities: epsilon=" << epsilon << " delta=" << delta << endl;
  }
  
  vector<vector<double> > RS07::chopZoneLikelihoods (int maxLen) const {
    const MID_HMM hmm (epsilon + (1-epsilon)*(1-2*delta), (1-epsilon)*delta, (1-epsilon)*delta,
		       (1-epsilon)*(1-2*delta), epsilon + (1-epsilon)*delta, (1-epsilon)*delta,
		       (1-epsilon)*(1-2*delta), epsilon + (1-epsilon)*delta, (1-epsilon)*delta,
		       verbose);
    return hmm.chopZoneLikelihoods (maxLen);
  }
}
