#pragma once

#include <vector>

using namespace std;

vector<double> implicit(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0,
    const double &step, const int &N, const double &eps, const unsigned &max_iter);

vector<double> trapeze(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0,
    const double &step, const int &N, const double &eps, const unsigned &max_iter);