#pragma once

#include <vector>

using namespace std;

vector<double> adams_bashforth_3(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0,
    const double &y1, const double &y2, const double &step, const int &N);

vector<double> adams_bashforth_4(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0,
    const double &y1, const double &y2, const double &y3, const double &step, const int &N);