#ifndef RK4_INCLUDED
#define RK4_INCLUDED

#include <functional>
#include <vector>

using namespace std;

typedef function<void(double,const vector<double>&,vector<double>&)> DerivFunc;   // params to DerivFunc are t,y,dy/dt

// A version of Runge-Kutta that uses a discretized time parameter
// Both the integer-valued time and the (computed) real time are passed to the derivative evaluator, which can use either
struct RungeKutta4 {
  const size_t dim;
  // scratch for Runge-Kutta
  vector<double> f0, f1, f2, f3, u1, u2, u3;
  RungeKutta4 (size_t);
  void step (DerivFunc evalDerivs, double t, double dt, const vector<double>& y0, vector<double>& y1);
};

#endif /* RK4_INCLUDED */
