#include "solvers.h"

vector<double> simple_explicit_euler(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, const double &step, const int &N)
{
   vector<double> y(N + 1);
   y[0] = y0;

   for (size_t i = 1; i < y.size(); i++)
      y[i] = y[i - 1] + step * f(grid[i], y[i - 1]);

   return y;
}

vector<double> modified_euler(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, const double &step, const int &N)
{
   // your implementation here...
   return vector<double>();
}

vector<double> improved_euler(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, const double &step, const int &N)
{
   // your implementation here...
   return vector<double>();
}
