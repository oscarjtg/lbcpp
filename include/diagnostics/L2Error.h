#ifndef L2ERRORDEF
#define L2ERRORDEF

#include "macroscopic/MacroscopicVariable.h"
#include "lattice/AbstractLattice.h"

#include <cmath>
#include <cassert>

/**
 * @brief Computes the L2 error norm for two scalar fields.
 */
template <typename T>
double ComputeL2ErrorScalar(MacroscopicVariable<T> &approx, MacroscopicVariable<T> &exact)
{
    assert(approx.GetNX() == exact.GetNX());
    assert(approx.GetNY() == exact.GetNY());
    assert(approx.GetNZ() == exact.GetNZ());

    double numerator = 0.0, denominator = 0.0;
    for (int k = 0; k < exact.GetNZ(); ++k)
    {
        for (int j = 0; j < exact.GetNY(); ++j)
        {
            for (int i = 0; i < exact.GetNX(); ++i)
            {
                double val_a = approx.GetValue(i, j, k);
                double val_e = exact.GetValue(i, j, k);
                double difference = val_a - val_e;
                numerator += difference * difference;
                denominator += val_e * val_e;
            }
        }
    }
    assert(denominator != 0);
    double l2_error = sqrt(numerator / denominator);
    return l2_error;
}

/**
 * @brief Computes the L2 error norm for a scalar field from an analytic function.
 */
template <typename T, typename Func>
double ComputeL2ErrorScalar(MacroscopicVariable<T> &array, Func analytic_function)
{
    double numerator = 0.0, denominator = 0.0;
    for (int k = 0; k < array.GetNZ(); ++k)
    {
        for (int j = 0; j < array.GetNY(); ++j)
        {
            for (int i = 0; i < array.GetNX(); ++i)
            {
                double val_a = array.GetValue(i, j, k);
                double val_e = analytic_function(i, j, k);
                double difference = val_a - val_e;
                numerator += difference * difference;
                denominator += val_e * val_e;
            }
        }
    }
    assert(denominator != 0);
    double l2_error = sqrt(numerator / denominator);
    return l2_error;
}

/**
 * @brief Computes the L2 error norm for two vector fields (ax, ay, az) and (ex, ey, ez).
 */
template <typename T>
double ComputeL2ErrorVector3(MacroscopicVariable<T> &ax,
                             MacroscopicVariable<T> &ay,
                             MacroscopicVariable<T> &az,
                             MacroscopicVariable<T> &ex,
                             MacroscopicVariable<T> &ey,
                             MacroscopicVariable<T> &ez)
{
    // Trust that the dimensions are the same.
    double numerator = 0.0, denominator = 0.0;
    for (int k = 0; k < ax.GetNZ(); ++k)
    {
        for (int j = 0; j < ax.GetNY(); ++j)
        {
            for (int i = 0; i < ax.GetNX(); ++i)
            {
                double val_ax = ax.GetValue(i, j, k);
                double val_ay = ay.GetValue(i, j, k);
                double val_az = az.GetValue(i, j, k);
                double val_ex = ex.GetValue(i, j, k);
                double val_ey = ey.GetValue(i, j, k);
                double val_ez = ez.GetValue(i, j, k);
                double diff_x = val_ax - val_ex;
                double diff_y = val_ay - val_ey;
                double diff_z = val_az - val_ez;
                numerator += diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
                denominator += val_ex * val_ex + val_ey * val_ey + val_ez * val_ez;
            }
        }
    }
    assert(denominator != 0);
    double l2_error = sqrt(numerator / denominator);
    return l2_error;
}

/**
 * @brief Computes the L2 error between two successive time steps of a set of Distribution Functions.
 */
template <typename T, int ND, int NQ>
double ComputeL2ErrorDistribution(AbstractLattice<T, ND, NQ> &f)
{
    double numerator = 0.0, denominator = 0.0;
    for (int k = 0; k < f.GetNZ(); ++k)
    {
        for (int j = 0; j < f.GetNY(); ++j)
        {
            for (int i = 0; i < f.GetNX(); ++i)
            {
                for (int q = 0; q < NQ; ++q)
                {
                    double f_curr = f.GetCurrFStar(q, i, j, k);
                    double f_prev = f.GetPrevFStar(q, i, j, k);
                    double diff = f_curr - f_prev;
                    numerator += diff * diff;
                    denominator += f_prev * f_prev;
                }
            }
        }
    }
    assert(denominator != 0);
    double l2_error = sqrt(numerator / denominator);
    return l2_error;
}

#endif // L2ERRORDEF