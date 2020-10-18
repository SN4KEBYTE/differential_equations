#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <direct.h>
#include "../utils/analytical.h"
#include "../utils/grid.h"
#include "../utils/excel.h"
#include "../euler_explicit/solvers.h"
#include "../euler_implicit/solvers.h"
#include "../runge_kutta/solvers.h"

using namespace std;

// task conditions
const double A = 0., B = 2., Y0 = 1.;

// for implicit methods
const double EPS = 1e-14;
const double MAX_ITER = 1000;

// steps
const vector<double> simple_steps = { 0.1, 0.05, 0.025 };
const vector<double> modified_steps = { 0.2, 0.1, 0.05 };
const vector<double> rk_steps = { 0.5, 0.25, 0.125, 0.1 };

// number of tests
const size_t SIMPLE_TEST_NUM = simple_steps.size();
const size_t MOD_TEST_NUM = modified_steps.size();
const size_t RK_TEST_NUM = rk_steps.size();

// root directory for resulting tables 
const string RES_DIR = "result/";

// directories for each type of methods
const string SIMPLE_DIR = RES_DIR + "simple/";
const string MOD_DIR = RES_DIR + "modified/";
const string RK_DIR = RES_DIR + "rk/";

// function
double f(const double &t, const double &y)
{
    return -25 * y + cos(t) + 25 * sin(t);
}

// analytical solution
double solution(const double &t)
{
    return exp(-25 * t) + sin(t);
}

double f2(const double &t, const double &y)
{
    return 2 * t * y;
}

double df2(const double &t, const double &y)
{
    return 2 * y;
}

// analytical solution
double solution2(const double &t)
{
    return exp(t * t);
}

int main()
{
    cout << "hello, world!" << endl;

    int N = (1. - 0.) / 0.05;
    auto grid = build_grid(0., 0.05, N);
    auto as = analytical_solution(grid, solution2);

    auto ns = implicit(grid, f2, df2, 1., 0.05, N, 1e-14, 1000);

    for (size_t i = 0; i < ns.size(); i++)
    {
        cout << as[i] << " " << ns[i] << " " << abs(as[i] - ns[i]) << endl;
    }

    system("pause");
    
    //_mkdir(RES_DIR.c_str());
   
    //_mkdir(SIMPLE_DIR.c_str());
    //_mkdir((SIMPLE_DIR + "simple/").c_str());
    //_mkdir((SIMPLE_DIR + "implicit/").c_str());
   
    //_mkdir(MOD_DIR.c_str());
    //_mkdir((MOD_DIR + "modified/").c_str());
    //_mkdir((MOD_DIR + "trapeze/").c_str());
    //
    //_mkdir(RK_DIR.c_str());

    //// SIMPLE
    //for (size_t i = 0; i < SIMPLE_TEST_NUM; i++)
    //{
    //    int N = (B - A) / simple_steps[i];
    //    auto grid = build_grid(A, simple_steps[i], N);
    //    auto as = analytical_solution(grid, solution);

    //    auto y_simple = simple(grid, f, Y0, simple_steps[i], N);
    //    auto y_implicit = implicit(grid, f, Y0, simple_steps[i], N, EPS, MAX_ITER);

    //    auto end = "h" + to_string(simple_steps[i]) + ".csv";

    //    ofstream out(SIMPLE_DIR + "simple/" + end);
    //    out.imbue(locale(""));
    //    dump_table(out, grid, y_simple, as);
    //    out.close();

    //    out.open(SIMPLE_DIR + "implicit/" + end);
    //    dump_table(out, grid, y_implicit, as);
    //    out.close();
    //}
    //
    //// MODIFIED
    //for (size_t i = 0; i < MOD_TEST_NUM; i++)
    //{
    //    int N = (B - A) / modified_steps[i];
    //    auto grid = build_grid(A, modified_steps[i], N);
    //    auto as = analytical_solution(grid, solution);

    //    auto y_mod = modified(grid, f, Y0, modified_steps[i], N);
    //    auto y_trapeze = trapeze(grid, f, Y0, modified_steps[i], N, EPS, MAX_ITER);

    //    auto end = "h" + to_string(modified_steps[i]) + ".csv";

    //    ofstream out(MOD_DIR + "modified/" + end);
    //    out.imbue(locale(""));
    //    dump_table(out, grid, y_mod, as);
    //    out.close();

    //    out.open(MOD_DIR + "trapeze/" + end);
    //    dump_table(out, grid, y_trapeze, as);
    //    out.close();
    //}

    //// RUNGE - KUTTA
    //for (size_t i = 0; i < RK_TEST_NUM; i++)
    //{
    //    int N = (B - A) / rk_steps[i];
    //    auto grid = build_grid(A, rk_steps[i], N);
    //    auto as = analytical_solution(grid, solution);

    //    auto y_rk = four_stage_runge_kutta(grid, f, Y0, rk_steps[i], N);

    //    auto end = "h" + to_string(rk_steps[i]) + ".csv";

    //    ofstream out(RK_DIR + end);
    //    out.imbue(locale(""));
    //    dump_table(out, grid, y_rk, as);
    //    out.close();
    //}

    return 0;
}