// solvers.cpp
#include "solvers.h"

vector<double> implicit(const vector<double> &grid, double(&f)(const double &, const double &), double(&df)(const double &, const double &),
    const double &y0, const double &step, const int &N, const double &eps, const unsigned &max_iter)
{
    auto y0_sol = modified(grid, f, y0, step, N);
    vector<double> y(N + 1);
    y[0] = y0;

    for (size_t i = 1; i < y.size(); i++)
    {
        double y_prev = y0_sol[i];
        double y_cur = y[i - 1] - (y[i - 1] + step * f(grid[i], y_prev) - y_prev) / (step * df(grid[i], y_prev) - 1);
        double accuracy = abs(y_cur - y_prev);
        
        unsigned iter = 0;

        while (accuracy > eps && iter <= max_iter)
        {
            y_prev = y_cur;
            y_cur = y[i - 1] - (y[i - 1] + step * f(grid[i], y_prev) - y_prev) / (step * df(grid[i], y_prev) - 1);
        
            accuracy = abs(y_cur - y_prev);
            iter++;
        }

        y[i] = y_cur;
    }

    return y;
}

vector<double> trapeze(const vector<double> &grid, double(&f)(const double &, const double &), double(&df)(const double &, const double &),
    const double &y0, const double &step, const int &N, const double &eps, const unsigned &max_iter)
{
    auto y0_sol = improved(grid, f, y0, step, N);
    vector<double> y(N + 1);
    y[0] = y0;

    for (size_t i = 1; i < y.size(); i++)
    {
        double y_prev = y0_sol[i];
        double y_cur = y[i - 1] - (y[i - 1] + (step / 2) * (f(grid[i - 1], y[i - 1]) + f(grid[i], y_prev)) - y_prev) /
            ((step / 2) * df(grid[i], y_prev) - 1);
        double accuracy = abs(y_cur - y_prev);
        
        unsigned iter = 0;

        while (accuracy > eps && iter <= max_iter)
        {
            y_prev = y_cur;
            y_cur = y[i - 1] - (y[i - 1] + (step / 2) * (f(grid[i - 1], y[i - 1]) + f(grid[i], y_prev)) - y_prev) /
                ((step / 2) * df(grid[i], y_prev) - 1);

            accuracy = abs(y_cur - y_prev);
            iter++;
        }

        y[i] = y_cur;
    }

    return y;
}