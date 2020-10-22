#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <direct.h>
#include "../utils/analytical.h"
#include "../utils/grid.h"
#include "../utils/excel.h"
#include "../adams_bashforth/solvers.h"
#include "../adams_moulton/solvers.h"
#include "../euler_explicit/solvers.h"
#include "../euler_implicit/solvers.h"
#include "../prediction_correction/solvers.h"
#include "../runge_kutta/solvers.h"

using namespace std;

#pragma region global const

// for implicit methods
const double EPS = 1e-14;
const double MAX_ITER = 1000;

// root directory for resulting tables 
const string RES_DIR = "result/";

#pragma endregion

#pragma region initial research data

// task conditions
const double I_A = 0., I_B = 1., I_Y0 = 1.;

// steps
const vector<double> I_STEPS = { 0.1, 0.05, 0.025 };

// root directory
const string INITIAL_DIR = "initial_result/";

// directories for each type of methods
const string I_AB3_DIR = INITIAL_DIR + "adams_bashforth3/";
const string I_AB4_DIR = INITIAL_DIR + "adams_bashforth4/";
const string I_AM3_DIR = INITIAL_DIR + "adams_moulton3/";
const string I_AM4_DIR = INITIAL_DIR + "adams_moulton4/";
const string I_EU_EXP_DIR = INITIAL_DIR + "euler_explicit/";
const string I_EU_IMP_DIR = INITIAL_DIR + "euler_implicit/";
const string I_PC_DIR = INITIAL_DIR + "prediction_correction/";
const string I_RK_DIR = INITIAL_DIR + "runge_kutta/";

// function 
double i_f(const double &t, const double &y)
{
    return 2 * t * y;
}

// function first derivative by y
double i_df(const double &t, const double &y)
{
    return 2 * t;
}

// analytical solution
double i_solution(const double &t)
{
    return exp(t * t);
}

#pragma endregion

//// task conditions
//const double A = 0., B = 2., Y0 = 1.;
//
//// steps
//const vector<double> simple_steps = { 0.1, 0.05, 0.025 };
//const vector<double> modified_steps = { 0.2, 0.1, 0.05 };
//const vector<double> rk_steps = { 0.5, 0.25, 0.125, 0.1 };
//
//// number of tests
//const size_t SIMPLE_TEST_NUM = simple_steps.size();
//const size_t MOD_TEST_NUM = modified_steps.size();
//const size_t RK_TEST_NUM = rk_steps.size();
//
//// root directory for resulting tables 
//const string RES_DIR = "result/";
//
//// directories for each type of methods
//const string SIMPLE_DIR = RES_DIR + "simple/";
//const string MOD_DIR = RES_DIR + "modified/";
//const string RK_DIR = RES_DIR + "rk/";
//
//// function
//double f(const double &t, const double &y)
//{
//    return -25 * y + cos(t) + 25 * sin(t);
//}
//
//// analytical solution
//double solution(const double &t)
//{
//    return exp(-25 * t) + sin(t);
//}
//
//double f2(const double &t, const double &y)
//{
//    return 2 * t * y;
//}
//
//double df2(const double &t, const double &y)
//{
//    return 2 * y;
//}
//
//// analytical solution
//double solution2(const double &t)
//{
//    return exp(t * t);
//}

int main()
{
#pragma region initial research
    // create neccesary directories
    _mkdir(INITIAL_DIR.c_str());

    _mkdir(I_AB3_DIR.c_str());
    _mkdir(I_AB4_DIR.c_str());
    _mkdir(I_AM3_DIR.c_str());
    _mkdir(I_AM4_DIR.c_str());
    _mkdir(I_EU_EXP_DIR.c_str());
    _mkdir(I_EU_IMP_DIR.c_str());
    _mkdir(I_PC_DIR.c_str());
    _mkdir(I_RK_DIR.c_str());

    // create ofstream and set proper locale to write numbers with comma
    ofstream out;
    out.imbue(locale(""));

    for (size_t i = 0; i < I_STEPS.size(); i++)
    {
        unsigned N = (I_B - I_A) / I_STEPS[i];
        auto grid = build_grid(I_A, I_STEPS[i], N);
        auto as = analytical_solution(grid, i_solution);
        auto end = "h" + to_string(I_STEPS[i]) + ".csv";

        // Runge-Kutta
        auto rk = four_stage_runge_kutta(grid, i_f, I_Y0, I_STEPS[i], N);
        out.open(I_RK_DIR + end);
        dump_table(out, grid, rk, as);
        out.close();

        // Adams-Bashforth
        auto ab3 = adams_bashforth_3(grid, i_f, I_Y0, rk[1], rk[2], I_STEPS[i], N);
        out.open(I_AB3_DIR + end);
        dump_table(out, grid, ab3, as);
        out.close();

        auto ab4 = adams_bashforth_4(grid, i_f, I_Y0, rk[1], rk[2], rk[3], I_STEPS[i], N);
        out.open(I_AB4_DIR + end);
        dump_table(out, grid, ab4, as);
        out.close();

        // Adams-Moulton
        auto am3 = adams_moulton_3(grid, i_f, i_df, I_Y0, I_STEPS[i], N, EPS, MAX_ITER);
        out.open(I_AM3_DIR + end);
        dump_table(out, grid, am3, as);
        out.close();

        auto am4 = adams_moulton_4(grid, i_f, i_df, I_Y0, I_STEPS[i], N, EPS, MAX_ITER);
        out.open(I_AM4_DIR + end);
        dump_table(out, grid, am4, as);
        out.close();

        // prediction-correction
        auto pc = prediction_correction(grid, i_f, I_Y0, rk[1], rk[2], I_STEPS[i], N);
        out.open(I_PC_DIR + end);
        dump_table(out, grid, pc, as);
        out.close();
    }

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