// solvers.cpp
#include "solvers.h"

vector<double> adams_bashforth_3(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0,
    const double &y1, const double &y2, const double &step, const int &N)
{
    vector<double> y(N + 1);
    y[0] = y0;
    y[1] = y1;
    y[2] = y2;

    for (size_t i = 3; i < y.size(); i++)
        y[i] = y[i - 1] + (step / 12) * (23 * f(grid[i - 1], y[i - 1]) - 16 * f(grid[i - 2], y[i - 2]) + 
            5 * f(grid[i - 3], y[i - 3]));

    return y;
}

vector<double> adams_bashforth_4(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0,
    const double &y1, const double &y2, const double &y3, const double &step, const int &N)
{
    vector<double> y(N + 1);
    y[0] = y0;
    y[1] = y1;
    y[2] = y2;
    y[3] = y3;

    for (size_t i = 4; i < y.size(); i++)
        y[i] = y[i - 1] + (step / 24) * (55 * f(grid[i - 1], y[i - 1]) - 59 * f(grid[i - 2], y[i - 2]) + 
            37 * f(grid[i - 3], y[i - 3]) - 9 * f(grid[i - 4], y[i - 4]));

    return y;
}