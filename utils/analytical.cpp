#include "analytical.h"

vector<double> analytical_solution(const vector<double> &grid, double(&f)(const double &))
{
    vector<double> as(grid.size());

    for (size_t i = 0; i < as.size(); i++)
        as[i] = f(grid[i]);

    return as;
}