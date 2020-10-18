// solvers.cpp
#include "solvers.h"

vector<double> four_stage_runge_kutta(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, 
    const double &step, const int &N)
{
    vector<double> y(N + 1);
    y[0] = y0;

    for (size_t i = 1; i < y.size(); i++)
    {
        double kn_1 = f(grid[i - 1], y[i - 1]);
        double kn_2 = f(grid[i - 1] + step / 2, y[i - 1] + step * kn_1 / 2);
        double kn_3 = f(grid[i - 1] + step / 2, y[i - 1] + step * kn_2 / 2);
        double kn_4 = f(grid[i - 1] + step, y[i - 1] + step * kn_3);

        double kn = (kn_1 + 2 * (kn_2 + kn_3) + kn_4) / 6;

        y[i] = y[i - 1] + step * kn;
    }

    return y;
}