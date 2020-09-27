#include "excel.h"

void dump_table(ofstream &f, const vector<double> &grid, const vector<double> &ns, const vector<double> &as)
{
    f << "t;y_числ;y_аналит;|y_числ - y_аналит|\n";

    for (size_t i = 0; i < grid.size(); i++)
        f << grid[i] << ";" << ns[i] << ";" << as[i] << ";" << abs(ns[i] - as[i]) << "\n";

    f.close();
}