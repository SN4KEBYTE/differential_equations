// solvers.cpp
#include "solvers.h"

vector<double> adams_moulton_3(const vector<double> &grid, double(&f)(const double &, const double &), double(&df)(const double &, const double &),
    const double &y0, const double &step, const int &N, const double &eps, const unsigned &max_iter)
{
    vector<double> y(N + 1);
    auto y0_sol = modified(grid, f, y0, step, N);
    y[0] = y0_sol[0];
    y[1] = y0_sol[1];

    for (size_t i = 2; i < y.size(); i++)
    {
        double y_prev = y0_sol[i];
        double y_cur = y_prev - (y[i - 1] + (step / 12) * (5 * f(grid[i], y_prev) + 8 * f(grid[i - 1], y[i - 1]) - f(grid[i - 2], y[i - 2])) - y_prev) /
            ((step / 12) * (5 * df(grid[i], y_prev)) - 1);
        double accuracy = abs(y_prev - y_cur);

        unsigned iter = 1;

        while (accuracy > eps && iter <= max_iter)
        {
            y_prev = y_cur;
            y_cur = y_prev - (y[i - 1] + (step / 12) * (5 * f(grid[i], y_prev) + 8 * f(grid[i - 1], y[i - 1]) - f(grid[i - 2], y[i - 2])) - y_prev) /
                ((step / 12) * (5 * df(grid[i], y_prev)) - 1);
            
            accuracy = abs(y_prev - y_cur);
            iter++;
        }

        y[i] = y_cur;
    }

    return y;
}

vector<double> adams_moulton_4(const vector<double> &grid, double(&f)(const double &, const double &), double(&df)(const double &, const double &),
    const double &y0, const double &step, const int &N, const double &eps, const unsigned &max_iter)
{
    vector<double> y(N + 1);
    auto y0_sol = four_stage_runge_kutta(grid, f, y0, step, N);
    y[0] = y0_sol[0];
    y[1] = y0_sol[1];
    y[2] = y0_sol[2];

    for (size_t i = 3; i < y.size(); i++)
    {
        double y_prev = y0_sol[i];
        double y_cur = y_prev - (y[i - 1] + (step / 24) * (9 * f(grid[i], y_prev) + 19 * f(grid[i - 1], y[i - 1]) - 5 * f(grid[i - 2], y[i - 2]) + f(grid[i - 3], y[i - 3])) - y_prev) /
            ((step / 24) * (9 * df(grid[i], y_prev)) - 1);
        double accuracy = abs(y_prev - y_cur);

        unsigned iter = 1;

        while (accuracy > eps && iter <= max_iter)
        {
            y_prev = y_cur;
            y_cur = y_prev - (y[i - 1] + (step / 24) * (9 * f(grid[i], y_prev) + 19 * f(grid[i - 1], y[i - 1]) - 5 * f(grid[i - 2], y[i - 2]) + f(grid[i - 3], y[i - 3])) - y_prev) /
                ((step / 24) * (9 * df(grid[i], y_prev)) - 1);

            accuracy = abs(y_prev - y_cur);
            iter++;
        }

        y[i] = y_cur;
    }

    return y;
}