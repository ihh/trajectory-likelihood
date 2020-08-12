#include <iostream>
#include "mid.h"

namespace TrajectoryLikelihood {

  MID_HMM::MID_HMM (double m2m, double m2i, double m2d, double i2m, double i2i, double i2d, double d2m, double d2i, double d2d, int verbose) :
    a(m2m), b(m2i), c(m2d),
    f(i2m), g(i2i), h(i2d),
    p(d2m), q(d2i), r(d2d)
  {
    if (verbose > 2)
      cerr << "Pair HMM probabilities: m2m=" << m2m << " m2i=" << m2i << " m2d=" << m2d << " i2m=" << i2m << " i2i=" << i2i << " d2m=" << d2m << " d2d=" << d2d << " d2i=" << d2i << endl;
    if (verbose > 3)
      cerr << "Pair HMM probability sums: m2*=" << (m2m+m2i+m2d) << " i2*=" << (i2m+i2i) << " d2*=" << (d2m+d2i+d2d) << endl;
  }

  vector<vector<double> > MID_HMM::chopZoneLikelihoods (int maxLen) const {
    vector<vector<double> > pGap (maxLen + 1, vector<double> (maxLen + 1, 0));
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
