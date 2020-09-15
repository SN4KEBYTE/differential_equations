#include <iostream>
#include <vector>
#include "grid.h"
#include "solvers.h"

using namespace std;

const double A = 0., B = 1., Y0 = 1.;

double f(const double &t, const double &y)
{
   return 2 * t * y;
}

int main()
{
   double step = 0.1;

   int N = (B - A) / step;
   auto g = build_grid(A, step, N);

   auto solution = simple_explicit_euler(g, f, 1.0, step, N);

   for (const auto &el : solution)
      cout << el << endl;

   cout << endl;

   system("pause");

   return 0;
}