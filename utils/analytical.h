// analytical.h
#pragma once

#include <vector>

using namespace std;

vector<double> analytical_solution(const vector<double> &grid, double(&f)(const double &));