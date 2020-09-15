#pragma once

#include <vector>

using namespace std;

vector<double> simple_explicit_euler(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, const double &step, const int &N);

vector<double> modified_euler(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, const double &step, const int &N);

vector<double> improved_euler(const vector<double> &grid, double(&f)(const double &, const double &), const double &y0, const double &step, const int &N);
