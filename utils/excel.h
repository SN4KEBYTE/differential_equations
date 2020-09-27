#pragma once

#include <fstream>
#include <vector>

using namespace std;

void dump_table(ofstream &f, const vector<double> &grid, const vector<double> &ns, const vector<double> &as);