#ifndef RK4_H
#define RK4_H
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;

class RK4
{
public:
    RK4();
};


class DPSolver
{
public:
    DPSolver();

    void InitValues(VectorXd const & initVals, double initTime)
    {
        m_values = initVals;
        m_time = initTime;
    }

    void CalcStep();

    double t() const {return m_time;}
    VectorXd vals() const {return m_values;}

protected:
    virtual VectorXd RecalcSystem(double time, VectorXd& val) = 0;

    VectorXd m_values, steps;
    double m_time;
    double m_step;
    double m_max;
    double m_min;
    Matrix<double, 7, 7> A;
    Matrix<double, 1, 7> b1, b2;
};

class MySolver : public DPSolver
{
public:
    double GoodSolution(double t) //решение корректное
    {
        double C = 0.5 + 0.3 / (2.1 * 2.1);
        return C * exp(-2.1 * t) + 0.3 / 2.1 * (t - 1 / 2.1);
    }
    VectorXd GoodSolution2(double t)
    {
        VectorXd ret(2);
        ret(0) = 2 * exp(-3 * t) - exp(-39 * t) + (1.0 / 3) * cos(t);
        ret(1) = -exp(-3 * t) + 2 * exp(-39 * t) - (1.0 / 3) * cos(t);
        return ret;
    }
protected:
    VectorXd RecalcSystem(double time, VectorXd& val) override
    {
        VectorXd ret(1);
        double y = val(0);
        ret(0) = 0.3 * time - 2.1 * y;
        return ret;
    }
};

class MySolver2 : public DPSolver
{
public:
    VectorXd GoodSolution2(double t)
    {
        VectorXd ret(2);
        ret(0) = 2 * exp(-3 * t) - exp(-39 * t) + (1.0 / 3) * cos(t);
        ret(1) = -exp(-3 * t) + 2 * exp(-39 * t) - (1.0 / 3) * cos(t);
        return ret;
    }
protected:
    VectorXd RecalcSystem(double time, VectorXd& val) override
    {
        VectorXd ret(2);
        double u1 = val(0);
        double u2 = val(1);
        ret(0) = 9 * u1 + 24 * u2 + 5 * cos(time) - (1.0 / 3) * sin(time);
        ret(1) = -24 * u1 - 51 * u2 - 9 * cos(time) + (1.0 / 3) * sin(time);
        return ret;
    }
};
#endif // RK4_H
