#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <locale>
#include <direct.h>
#include "../utils/grid.h"
#include "solvers.h"
#include "../utils/analytical.h"
#include "../utils/excel.h"

using namespace std;

// task conditions
const double A = 0., B = 1., Y0 = 1.;

// number of tests
const size_t TEST_NUM = 3;

// root directory for resulting tables 
const string RES_DIR = "result/";

// function
double f(const double &t, const double &y)
{
    return 2 * t * y;
}

double solution(const double &t)
{
    return exp(t * t);
}

int main()
{
    // initialize step sizes
    const vector<double> steps = { 0.1, 0.05, 0.025 };

    // directories for each method's resulting tables
    const string IMP_DIR = RES_DIR + "imp/";
    const string TR_DIR = RES_DIR + "tr/";

    // create directories for resulting tables
    _mkdir(RES_DIR.c_str());
    _mkdir(IMP_DIR.c_str());
    _mkdir(TR_DIR.c_str());

    for (size_t i = 0; i < TEST_NUM; i++)
    {
        // calculate number of segments, build grid and calculate analytical solution
        int N = (B - A) / steps[i];
        auto grid = build_grid(A, steps[i], N);
        auto as = analytical_solution(grid, solution);

        // calculate solution using all methods for current grid
        auto y_imp = implicit(grid, f, Y0, steps[i], N, steps[i] / 10, 1000);
        auto y_tr = trapeze(grid, f, Y0, steps[i], N, steps[i] / 10, 1000);

        // dump results to files
        auto end = "h" + to_string(steps[i]) + ".csv";

        ofstream out(IMP_DIR + end);
        out.imbue(locale(""));
        dump_table(out, grid, y_imp, as);
        out.close();

        out.open(TR_DIR + end);
        dump_table(out, grid, y_tr, as);
        out.close();
    }   

    std::cout << "DONE" << endl;
    std::system("pause");

    return 0;
}
