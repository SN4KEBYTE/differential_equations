// solvers.h
#pragma once

#include <vector>

using namespace std;

vector<double> simple(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, const double &step, const int &N);

vector<double> modified(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, const double &step, const int &N);

vector<double> improved(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, const double &step, const int &N);
