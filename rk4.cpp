#include <cmath>
#include "rk4.h"

RungeKutta4::RungeKutta4 (size_t dim) :
  dim (dim),
  f0 (dim), f1 (dim), f2 (dim), f3 (dim),
  u1 (dim), u2 (dim), u3 (dim)
{ }

void RungeKutta4::step (DerivFunc evalDerivs, double t, double dt, const vector<double>& y0, vector<double>& y1) {
  // adapted from John Burkardt's rk4: https://people.sc.fsu.edu/~jburkardt/cpp_src/rk4/rk4.cpp
  size_t i;
  const double t0 = t, t1 = t + dt/2, t2 = t + dt/2, t3 = t + dt;
  
  evalDerivs (t0, y0, f0);

  for (i = 0; i < dim; ++i)
    u1[i] = y0[i] + dt * f0[i] / 2;
  evalDerivs (t1, u1, f1);

  for (i = 0; i < dim; ++i)
    u2[i] = y0[i] + dt * f1[i] / 2;
  evalDerivs (t2, u2, f2);

  for (i = 0; i < dim; ++i)
    u3[i] = y0[i] + dt * f2[i];
  evalDerivs (t3, u3, f3);

  for (i = 0; i < dim; ++i) {
    y1[i] = y0[i] + dt * (f0[i] + 2*f1[i] + 2*f2[i] + f3[i]) / 6;
    if (isnan(y1[i]))
      throw runtime_error ("rk4 isNaN");
    if (isinf(y1[i]))
      throw runtime_error ("rk4 isInf");
  }
}
