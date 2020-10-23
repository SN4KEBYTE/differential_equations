// solvers.h
#pragma once

#include <vector>
#include "../runge_kutta/solvers.h"

using namespace std;

vector<double> adams_moulton_3(const vector<double> &grid, double(&f)(const double &, const double &), double(&df)(const double &, const double &),
    const double &y0, const double &step, const int &N, const double &eps, const unsigned &max_iter);

vector<double> adams_moulton_4(const vector<double> &grid, double(&f)(const double &, const double &), double(&df)(const double &, const double &),
    const double &y0, const double &step, const int &N, const double &eps, const unsigned &max_iter);