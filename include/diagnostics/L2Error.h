#ifndef L2ERRORDEF
#define L2ERRORDEF

#include "macroscopic/AllFields.h"
#include "lattice/AbstractLattice.h"

#include <cmath>
#include <cassert>

/**
 * @brief Computes the L2 error norm for two scalar fields.
 */
template <typename T>
double ComputeL2ErrorScalar(ScalarField<T> &approx, ScalarField<T> &exact)
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
double ComputeL2ErrorScalar(ScalarField<T> &array, Func analytic_function, double subtract=0)
{
    double numerator = 0.0, denominator = 0.0;
    for (int k = 0; k < array.GetNZ(); ++k)
    {
        for (int j = 0; j < array.GetNY(); ++j)
        {
            for (int i = 0; i < array.GetNX(); ++i)
            {
                double val_a = array.GetValue(i, j, k) - subtract;
                double val_e = analytic_function(i, j, k) - subtract;
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
 * @brief Computes the L2 error norm for a component (int a) of a VectorField, 
 * compared to an analytic_function(i, j, k) for that component.
 */
template <typename T, typename Func>
double ComputeL2ErrorVectorComponent(VectorField<T> &array, int a, Func analytic_function)
{
    double numerator = 0.0, denominator = 0.0;
    for (int k = 0; k < array.GetNZ(); ++k)
    {
        for (int j = 0; j < array.GetNY(); ++j)
        {
            for (int i = 0; i < array.GetNX(); ++i)
            {
                double val_a = array.GetValue(a, i, j, k);
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
 * @brief Computes the L2 error norm for a component (int a, int b) of a TensorField, 
 * compared to an analytic_function(i, j, k) for that component.
 */
template <typename T, typename Func>
double ComputeL2ErrorTensorComponent(TensorField<T> &array, int a, int b, Func analytic_function)
{
    double numerator = 0.0, denominator = 0.0;
    for (int k = 0; k < array.GetNZ(); ++k)
    {
        for (int j = 0; j < array.GetNY(); ++j)
        {
            for (int i = 0; i < array.GetNX(); ++i)
            {
                double val_a = array.GetValue(a, b, i, j, k);
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
 * @brief Computes the L2 error norm for two vector fields (ax, ay, az) and (ex, ey, ez),
 * given as ScalarField instances.
 */
template <typename T>
double ComputeL2ErrorVector3(ScalarField<T> &ax,
                             ScalarField<T> &ay,
                             ScalarField<T> &az,
                             ScalarField<T> &ex,
                             ScalarField<T> &ey,
                             ScalarField<T> &ez)
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
 * @brief Computes the L2 error norm for two VectorFields (ax, ay, az) and (ex, ey, ez),
 * 
 */
template <typename T>
double ComputeL2ErrorVector3(VectorField<T>& approx, VectorField<T>& exact)
{
    // Trust that the dimensions are the same.
    double numerator = 0.0, denominator = 0.0;
    for (int k = 0; k < approx.GetNZ(); ++k)
    {
        for (int j = 0; j < approx.GetNY(); ++j)
        {
            for (int i = 0; i < approx.GetNX(); ++i)
            {
                double val_ax = approx.GetValue(0, i, j, k);
                double val_ay = approx.GetValue(1, j, k);
                double val_az = approx.GetValue(2, j, k);
                double val_ex = exact.GetValue(0, i, j, k);
                double val_ey = exact.GetValue(1, i, j, k);
                double val_ez = exact.GetValue(2, i, j, k);
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
                    denominator += f_curr * f_curr;
                }
            }
        }
    }
    assert(denominator != 0);
    double l2_error = sqrt(numerator / denominator);
    return l2_error;
}

#endif // L2ERRORDEF