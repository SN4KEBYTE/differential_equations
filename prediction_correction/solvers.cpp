// solvers.cpp
#include "solvers.h"

vector<double> prediction_correction_3(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0,
    const double &y1, const double &y2, const double &step, const int &N)
{
    auto y = adams_bashforth_3(grid, f, y0, y1, y2, step, N);

    for (size_t i = 2; i < y.size(); i++)
        y[i] = y[i - 1] + (step / 12) * (5 * f(grid[i], y[i]) + 8 * f(grid[i - 1], y[i - 1]) - f(grid[i - 2], y[i - 2]));

    return y;
}

vector<double> prediction_correction_4(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0,
    const double &y1, const double &y2, const double &y3, const double &step, const int &N)
{
    auto y = adams_bashforth_4(grid, f, y0, y1, y2, y3, step, N);

    for (size_t i = 3; i < y.size(); i++)
        y[i] = y[i - 1] + (step / 24) * (9 * f(grid[i], y[i]) + 19 * f(grid[i - 1], y[i - 1]) - 5 * f(grid[i - 2], y[i - 2]) + f(grid[i - 3], y[i - 3]));

    return y;
}