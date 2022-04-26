    /* Code inspired by sobol.cc available at:
       https://web.maths.unsw.edu.au/~fkuo/sobol/
    */

#define _USE_MATH_DEFINES

#include <cmath>   
#include <iostream>
#include <iomanip>

using namespace std;

double sobol_2D_integral(unsigned N, double (*fptr)(double, double)) { 
  /* Estimate the 2 dimensional integral over the unit square of function provided as second
     argument. Uses the first N Sobol sequence points as the empirical distribution approximation
     to the uniform distribution over the unit square. Estimates the Expectation of f using
     this distribution.
  */
  // L = max number of bits needed for this many points:
  unsigned L = (unsigned)ceil(log((double)N) / log(2.0));

  // C[i] = index from the right of the first zero bit of i < N:
  unsigned* C = new unsigned[N];
  C[0] = 1;
  for (unsigned i = 1; i <= N - 1; i++) {
    C[i] = 1;
    unsigned value = i;
    while (value & 1) {
      value >>= 1;
      C[i]++;
    }
  }

  // Compute direction numbers V_x[1] to V_x[L], scaled by pow(2,32)
  unsigned* V_x = new unsigned[L + 1];
  for (unsigned i = 1; i <= L; i++)
    V_x[i] = 1 << (32 - i);

  // Compute direction numbers V_y[1] to V_y[L]:
  unsigned* V_y = new unsigned[L + 1];
  V_y[1] = 1 << 31;
  for (unsigned i = 2; i <= L; i++) {
    V_y[i] = V_y[i - 1] ^ (V_y[i - 1] >> 1);
  }

  // Use V_x, V_y and C to get (x, y) Sobol points.
  // Accumulate f(x, y) into accum for mean estimate:
  double accum{ fptr(0,0) }; // Initialize with first point of Sobol sequence!
  unsigned X{ 0 }, Y{ 0 };
  for (unsigned i = 1; i <= N - 1; i++) {
    X = X ^ V_x[C[i - 1]];
    Y = Y ^ V_y[C[i - 1]];
    double x = static_cast<double>(X) / pow(2.0, 32);
    double y = static_cast<double>(Y) / pow(2.0, 32);
    double val = fptr(x, y);
    accum += val;
  }

  double est_mean = accum / N;

  delete[] V_x;
  delete[] V_y;
  delete[] C;

  return est_mean;
}

inline double circle4Indicator(double x, double y) {
  return x*x + y*y < 1.0 ? 4.0 : 0.0;
}


int main() {
  long a[] = { 100, 1000, 10000, 100000, 1000000 };
  double est_pi{ 0 };

  cout << "MC simulate fraction of (x,y) 2-D Sobol sequence points falling inside unit circle:\n";
  for (long n : a) {
    //double Pi_calc = sobol_2D_integral(n, [](double x, double y) { return x*x + y*y < 1.0 ? 4.0 : 0.0;} );
    double Pi_calc = sobol_2D_integral(n, circle4Indicator);
    cout << setprecision(6) << fixed;
    cout << "N: " << setw(10) << n << ", Calculated Pi=" << Pi_calc
      << ", Diff=" << Pi_calc - M_PI << endl;
  }

}
