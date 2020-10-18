// solvers.cpp
#include "solvers.h"
#include "../euler_explicit/solvers.h"

vector<double> implicit(const vector<double> &grid, double(&f)(const double &, const double &), double(&df)(const double &, const double &), 
    const double &y0, const double &step, const int &N, const double &eps, const unsigned &max_iter)
{
    auto y0_sol = modified(grid, f, y0, step, N);
    vector<double> y(N + 1);
    y[0] = y0;

    for (size_t i = 1; i < y.size(); i++)
    {
        double y_prev = y0_sol[i];
        double y_cur = y[i - 1] - (y[i - 1] + step * f(grid[i], y_prev) - y_prev) / (step * df(grid[i], y_prev));
        unsigned iter = 0;

        do
        {
            y_prev = y_cur;
            y_cur = y[i - 1] - (y[i - 1] + step * f(grid[i], y_prev) - y_prev) / (step * df(grid[i], y_prev));
            iter++;
        } while (abs(y_cur - y_prev) > eps && iter <= max_iter);

        y[i] = y_cur;

        /*unsigned iter = 0;

        double y_prev = y[i - 1];
        double y_cur = y[i - 1] + step * f(grid[i], y_prev);

        do
        {
            y_prev = y_cur;
            y_cur = y[i - 1] + step * f(grid[i], y_prev);
            iter++;
        } while (abs(y_cur - y_prev) >= eps && iter <= max_iter);

        y[i] = newton_solver(y0_sol[i], grid[i], f, df, eps, max_iter);*/
    }

    return y;
}

vector<double> trapeze(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0,
    const double &step, const int &N, const double &eps, const unsigned &max_iter)
{
    auto y0_sol = modified(grid, f, y0, step, N);
    vector<double> y(N + 1);
    y[0] = y0;

    for (size_t i = 1; i < y.size(); i++)
    {
        /*unsigned iter = 0;

        double y_prev = y[i - 1];
        double y_cur = y[i - 1] + (step / 2) * (f(grid[i - 1], y[i - 1]) + f(grid[i], y_prev));

        do
        {
            y_prev = y_cur;
            y_cur = y[i - 1] + (step / 2) * (f(grid[i - 1], y[i - 1]) + f(grid[i], y_prev));
            iter++;
        } while (abs(y_cur - y_prev) >= eps && iter <= max_iter);

        y[i] = newton_solver(y0_sol[i], grid[i], f, df, eps, max_iter);*/
    }

    return y;
}