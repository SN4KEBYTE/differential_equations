// grid.cpp
#include "grid.h"

vector<double> build_grid(const double &a, const double &step, const int &N)
{
    vector<double> grid(N + 1);

    for (size_t i = 0; i < N + 1; i++)
        grid[i] = a + i * step;

    return grid;
}
