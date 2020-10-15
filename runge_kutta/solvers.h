#pragma once

#include <vector>

using namespace std;

vector<double> four_stage_runge_kutta(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, 
    const double &step, const int &N);