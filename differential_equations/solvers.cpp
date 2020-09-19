#include "solvers.h"

vector<double> simple_explicit_euler(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, const double &step, const int &N)
{
    vector<double> y(N + 1);
    y[0] = y0;

    for (size_t i = 1; i < y.size(); i++)
        y[i] = y[i - 1] + step * f(grid[i - 1], y[i - 1]);

    return y;
}

vector<double> modified_euler(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, const double &step, const int &N)
{
    vector<double> y(N + 1);
    y[0] = y0;

    for (size_t i = 1; i < y.size(); i++)
        y[i] = y[i - 1] + (step / 2) * (f(grid[i - 1], y[i - 1]) + f(grid[i], y[i - 1] + (step / 2) * f(grid[i - 1], y[i - 1])));

    return y;
}

vector<double> improved_euler(const vector<double> &grid, double(&f)(const double &, const double &), const double &z0, const double &step, const int &N)
{
    vector<double> z(N + 1);
    z[0] = z0;

    for (size_t i = 1; i < z.size(); i++)
	z[i] = z[i - 1] + step * f(grid[i+1/2], z[i - 1]+(step/2)*f(grid[i], z[i - 1]));

    return z;
}
