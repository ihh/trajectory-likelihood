#ifndef CIM_INCLUDED
#define CIM_INCLUDED

#include <vector>
#include "trajec.h"

namespace TrajectoryLikelihood {

  using namespace std;
  
  struct CumulativeIndelModel {
    const IndelParams& params;
    const double ri, rd, gi, gd;
    vector<vector<double> > A;  // Ai, Ad, Aid
    double dt, tMax;
    int verbose;
    CumulativeIndelModel (const IndelParams&, double tMax, double dt, int verbose);
    double Lit (double t) const;
    double Pmt (double t) const;
    double Pit (const vector<double>& A, double t) const;
    double Pdt (const vector<double>& A, double t) const;
    double Pidt (const vector<double>& A, double t) const;
    double git (const vector<double>& A, double t) const;
    double gdt (const vector<double>& A, double t) const;

    vector<vector<double> > chopZoneLikelihoods (int maxLen) const;
  };
}

#endif /* CIM_INCLUDED */
