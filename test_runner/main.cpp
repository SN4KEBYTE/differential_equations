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
const string I_PC3_DIR = INITIAL_DIR + "prediction_correction3/";
const string I_PC4_DIR = INITIAL_DIR + "prediction_correction4/";
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

#pragma region convergence research data

// task conditions
const double C_A = 0., C_B = 2., C_Y0 = 1.;

// steps
const vector<double> C_SIMPLE_STEPS = { 0.1, 0.05, 0.025 };
const vector<double> C_MODIFIED_STEPS = { 0.2, 0.1, 0.05 };
const vector<double> C_RK_STEPS = { 0.5, 0.25, 0.125, 0.1 };

// number of tests
const size_t SIMPLE_TEST_NUM = C_SIMPLE_STEPS.size();
const size_t MOD_TEST_NUM = C_MODIFIED_STEPS.size();
const size_t RK_TEST_NUM = C_RK_STEPS.size();

// root directory for resulting tables 
const string CONV_DIR = "convergence_result/";

// directories for each type of methods
const string SIMPLE_DIR = CONV_DIR + "simple/";
const string MOD_DIR = CONV_DIR + "modified/";
const string RK_DIR = CONV_DIR + "rk/";

// function
double c_f(const double &t, const double &y)
{
    return -25 * y + cos(t) + 25 * sin(t);
}

// function first derivative by y
double c_df(const double &t, const double &y)
{
    return -25.;
}

// analytical solution
double c_sol(const double &t)
{
    return exp(-25 * t) + sin(t);
}

#pragma endregion

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
    _mkdir(I_PC3_DIR.c_str());
    _mkdir(I_PC4_DIR.c_str());
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
        auto pc3 = prediction_correction_3(grid, i_f, I_Y0, rk[1], rk[2], I_STEPS[i], N);
        out.open(I_PC3_DIR + end);
        dump_table(out, grid, pc3, as);
        out.close();

        auto pc4 = prediction_correction_4(grid, i_f, I_Y0, rk[1], rk[2], rk[3], I_STEPS[i], N);
        out.open(I_PC4_DIR + end);
        dump_table(out, grid, pc4, as);
        out.close();

        // implicit Euler
        auto imp = implicit(grid, i_f, i_df, I_Y0, I_STEPS[i], N, EPS, MAX_ITER);
        out.open(I_EU_IMP_DIR + end);
        dump_table(out, grid, imp, as);
        out.close();
    }
#pragma endregion

#pragma region convergence research
    _mkdir(CONV_DIR.c_str());

    _mkdir(SIMPLE_DIR.c_str());
    _mkdir((SIMPLE_DIR + "explicit/").c_str());
    _mkdir((SIMPLE_DIR + "implicit/").c_str());

    _mkdir(MOD_DIR.c_str());
    _mkdir((MOD_DIR + "modified/").c_str());
    _mkdir((MOD_DIR + "trapeze/").c_str());
    
    _mkdir(RK_DIR.c_str());

    // SIMPLE
    for (size_t i = 0; i < SIMPLE_TEST_NUM; i++)
    {
        int N = (C_B - C_A) / C_SIMPLE_STEPS[i];
        auto grid = build_grid(C_A, C_SIMPLE_STEPS[i], N);
        auto as = analytical_solution(grid, c_sol);

        auto y_simple = simple(grid, c_f, C_Y0, C_SIMPLE_STEPS[i], N);
        auto y_implicit = implicit(grid, c_f, c_df, C_Y0, C_SIMPLE_STEPS[i], N, EPS, MAX_ITER);

        auto end = "h" + to_string(C_SIMPLE_STEPS[i]) + ".csv";

        ofstream out(SIMPLE_DIR + "explicit/" + end);
        out.imbue(locale(""));
        dump_table(out, grid, y_simple, as);
        out.close();

        out.open(SIMPLE_DIR + "implicit/" + end);
        dump_table(out, grid, y_implicit, as);
        out.close();
    }
    
    // MODIFIED
    for (size_t i = 0; i < MOD_TEST_NUM; i++)
    {
        int N = (C_B - C_A) / C_MODIFIED_STEPS[i];
        auto grid = build_grid(C_A, C_MODIFIED_STEPS[i], N);
        auto as = analytical_solution(grid, c_sol);

        auto y_mod = modified(grid, c_f, C_Y0, C_MODIFIED_STEPS[i], N);
        auto y_trapeze = trapeze(grid, c_f, c_df, C_Y0, C_MODIFIED_STEPS[i], N, EPS, MAX_ITER);

        auto end = "h" + to_string(C_MODIFIED_STEPS[i]) + ".csv";

        ofstream out(MOD_DIR + "modified/" + end);
        out.imbue(locale(""));
        dump_table(out, grid, y_mod, as);
        out.close();

        out.open(MOD_DIR + "trapeze/" + end);
        dump_table(out, grid, y_trapeze, as);
        out.close();
    }

    // RUNGE - KUTTA
    for (size_t i = 0; i < RK_TEST_NUM; i++)
    {
        int N = (C_B - C_A) / C_RK_STEPS[i];
        auto grid = build_grid(C_A, C_RK_STEPS[i], N);
        auto as = analytical_solution(grid, c_sol);

        auto y_rk = four_stage_runge_kutta(grid, c_f, C_Y0, C_RK_STEPS[i], N);

        auto end = "h" + to_string(C_RK_STEPS[i]) + ".csv";

        ofstream out(RK_DIR + end);
        out.imbue(locale(""));
        dump_table(out, grid, y_rk, as);
        out.close();
    }
#pragma endregion

    return 0;
}