#include "rk4.h"

using namespace std;


RK4::RK4()
{
    cout << "all good";
}


DPSolver::DPSolver()
{
    // Инициализация коэффициентов метода Рунге-Кутта
    //Коэффициенты из метода Дорманда-Принца
    steps.resize(7);
    A << 0, 0, 0, 0, 0, 0, 0,
            1.0/5, 0, 0, 0, 0, 0, 0,
            3.0/40, 9.0/40, 0, 0, 0, 0, 0,
            44.0/45, -56.0/15, 32.0/9, 0, 0, 0, 0,
            19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729, 0, 0, 0,
            9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656, 0, 0,
            35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84, 0;

    b1 << 35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84, 0;
    b2 << 5179.0/57600, 0, 7571.0/16695, 393.0/640, -92097.0/339200, 187.0/2100, 1.0/40;

    steps << 0, 1.0/5, 3.0/10, 4.0/5, 8.0/9, 1, 1;

    m_step = 0.01;
    m_max = 0.0001;
    m_min = 0.0000001;
}

void DPSolver::CalcStep()
{
    MatrixXd tmp(7, m_values.size());
    VectorXd x1, x2;
    double diff = 1;
    do
    {
        VectorXd val = m_values;
        tmp.row(0) = RecalcSystem(m_time, val);
        for (int i = 1; i < 7; i++)
        {
            VectorXd val = m_values;
            for (int j = 0; j < i; j++)
                val += tmp.row(j) * A(i,j) * m_step;
            tmp.row(i) = RecalcSystem(m_time + steps(i) * m_step, val);
        }
        x1 = b1 * tmp;
        x2 = b2 * tmp;
        diff = (x1 - x2).cwiseAbs().maxCoeff() * m_step;
        if (diff > m_max)
            m_step /= 2;
        if (diff < m_min)
            m_step *= 2;
        std::cout << "Шаг: " << m_step << endl;
    }

    while (diff > m_max);
    m_values = m_values + x1 * m_step;
    m_time += m_step;
}
