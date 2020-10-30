// solvers.h
#pragma once

#include <vector>
#include "../adams_bashforth/solvers.h"

using namespace std;

vector<double> prediction_correction_3(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0,
    const double &y1, const double &y2, const double &step, const int &N);

vector<double> prediction_correction_4(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0,
    const double &y1, const double &y2, const double &y3, const double &step, const int &N);