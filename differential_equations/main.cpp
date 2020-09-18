#include <iostream>
#include <vector>
#include <locale>
#include "grid.h"
#include "solvers.h"

using namespace std;

const double A = 0., B = 1., Y0 = 1.;

double f(const double &t, const double &y)
{
    return 2 * t * y;
}

int main()
{
    cout.imbue(locale(""));

    double step = 0.025;

    int N = (B - A) / step;
    auto g = build_grid(A, step, N);

    cout << "Total elements: " << N << endl;

    for (const auto &el : g)
        cout << el << endl;

    /*auto solution = simple_explicit_euler(g, f, 1.0, step, N);

    for (const auto &el : solution)
       cout << el << endl;

    cout << endl;*/
    
    auto solution_3 = improved_euler(g, f, 1.0, step, N);
    
    for (const auto& el : solution_3)
	cout << el << endl;
    cout << endl;

    system("pause");

    return 0;
}
