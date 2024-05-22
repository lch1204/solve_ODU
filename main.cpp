#include <QCoreApplication>
#include "rk4.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

int main()
{
    VectorXd initVals(1);
    initVals(0) = 0.5; ///начальное значение y(0)=d=0.5

    MySolver solve;
    solve.InitValues(initVals, 0);

    double maxDiscrepancy1 = 0.0;
    double maxDiscrepancy2 = 0.0;

    while (solve.t() < 5) {
            solve.CalcStep();
            double t = solve.t();
            double numericalSolution = solve.vals()(0);
            double goodSolution = solve.GoodSolution(t);
            double discrepancy = std::abs(goodSolution - numericalSolution);
            maxDiscrepancy1 = std::max(maxDiscrepancy1, discrepancy);

            cout << "t: " << t << ", y(t): " << numericalSolution << ", u(t): " << goodSolution << ", |u(t) - y(t)|: " << discrepancy << endl;
    }



    VectorXd initVals2(2);
    initVals2(0) = 4/3;
    initVals2(1) = 2/3;

    MySolver2 solve2;
    solve2.InitValues(initVals2, 0);



    while (solve2.t() < 5) {
            solve2.CalcStep();
            double t = solve2.t();
            VectorXd numericalSolution = solve2.vals();
            VectorXd goodSolution2 = solve2.GoodSolution2(t);
            double discrepancy = (numericalSolution - goodSolution2).cwiseAbs().maxCoeff();
            maxDiscrepancy2 = std::max(maxDiscrepancy2, discrepancy);

            cout << "t: " << t << ", u1(t): " << numericalSolution(0) << ", y1(t): " << goodSolution2(0) <<", u2(t): " << numericalSolution(1)
                     << ", y2(t): " << goodSolution2(1) << ", разность: " << discrepancy << endl;
        }
        cout << "Максимальная разность для первого уравнения: " << maxDiscrepancy1 << endl;
        cout << "Максимальная разность системы уравнений: " << maxDiscrepancy2 << endl;

    return 0;
}
