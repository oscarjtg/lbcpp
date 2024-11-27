#ifndef FINITEDIFFERENCEDEF
#define FINITEDIFFERENCEDEF

#include "macroscopic/MacroscopicVariable.h"

#include <cmath>

template <typename T>
double ComputeCentredDiff(const ScalarField<T>& M, int i, int j, int k, int dx, int dy, int dz)
{
    // Distance between centre and forward point.
    double h = sqrt(static_cast<double>(dx*dx + dy*dy + dz*dz));

    // Forward point.
    double uf = M.GetValueWrap(i+dx, j+dy, k+dz);

    // Backward point.
    double ub = M.GetValueWrap(i-dx, j-dy, k-dz);

    // Centred difference approximation to first derivative.
    return (uf - ub) / (2.0*h);
}

template <typename T>
double ComputeForwardDiff(const ScalarField<T>& M, int i, int j, int k, int dx, int dy, int dz)
{
    // Distance between centre and forward point.
    double h = sqrt(static_cast<double>(dx*dx + dy*dy + dz*dz));

    // Forward point.
    double uf = M.GetValueWrap(i+dx, j+dy, k+dz);

    // Centre point.
    double uc = M.GetValue(i, j, k);

    // Forward different approximation to first derivative.
    return (uf - uc) / h;
}

template <typename T>
double ComputeBackwardDiff(const ScalarField<T>& M, int i, int j, int k, int dx, int dy, int dz)
{
    // Distance between centre and forward point.
    double h = sqrt(static_cast<double>(dx*dx + dy*dy + dz*dz));

    // Backwards point.
    double ub = M.GetValueWrap(i-dx, j-dy, k-dz);

    // Centre point.
    double uc = M.GetValue(i, j, k);

    // Forward different approximation to first derivative.
    return (uc - ub) / h;
}

template <typename T>
double ComputeCentredDiff(const VectorField<T>& M, int a, int i, int j, int k, int dx, int dy, int dz)
{
    // Distance between centre and forward point.
    double h = sqrt(static_cast<double>(dx*dx + dy*dy + dz*dz));

    // Forward point.
    double uf = M.GetValueWrap(a, i+dx, j+dy, k+dz);

    // Backward point.
    double ub = M.GetValueWrap(a, i-dx, j-dy, k-dz);

    // Centred difference approximation to first derivative.
    return (uf - ub) / (2.0*h);
}

template <typename T>
double ComputeForwardDiff(const VectorField<T>& M, int a, int i, int j, int k, int dx, int dy, int dz)
{
    // Distance between centre and forward point.
    double h = sqrt(static_cast<double>(dx*dx + dy*dy + dz*dz));

    // Forward point.
    double uf = M.GetValueWrap(a, i+dx, j+dy, k+dz);

    // Centre point.
    double uc = M.GetValue(a, i, j, k);

    // Forward different approximation to first derivative.
    return (uf - uc) / h;
}

template <typename T>
double ComputeBackwardDiff(const VectorField<T>& M, int a, int i, int j, int k, int dx, int dy, int dz)
{
    // Distance between centre and forward point.
    double h = sqrt(static_cast<double>(dx*dx + dy*dy + dz*dz));

    // Backwards point.
    double ub = M.GetValueWrap(a, i-dx, j-dy, k-dz);

    // Centre point.
    double uc = M.GetValue(a, i, j, k);

    // Forward different approximation to first derivative.
    return (uc - ub) / h;
}

template <typename T>
double ComputeCentredDiff(const TensorField<T>& M, int a, int b, int i, int j, int k, int dx, int dy, int dz)
{
    // Distance between centre and forward point.
    double h = sqrt(static_cast<double>(dx*dx + dy*dy + dz*dz));

    // Forward point.
    double uf = M.GetValueWrap(a, b, i+dx, j+dy, k+dz);

    // Backward point.
    double ub = M.GetValueWrap(a, b, i-dx, j-dy, k-dz);

    // Centred difference approximation to first derivative.
    return (uf - ub) / (2.0*h);
}

template <typename T>
double ComputeForwardDiff(const TensorField<T>& M, int a, int b, int i, int j, int k, int dx, int dy, int dz)
{
    // Distance between centre and forward point.
    double h = sqrt(static_cast<double>(dx*dx + dy*dy + dz*dz));

    // Forward point.
    double uf = M.GetValueWrap(a, b, i+dx, j+dy, k+dz);

    // Centre point.
    double uc = M.GetValue(a, b, i, j, k);

    // Forward different approximation to first derivative.
    return (uf - uc) / h;
}

template <typename T>
double ComputeBackwardDiff(const TensorField<T>& M, int a, int b, int i, int j, int k, int dx, int dy, int dz)
{
    // Distance between centre and forward point.
    double h = sqrt(static_cast<double>(dx*dx + dy*dy + dz*dz));

    // Backwards point.
    double ub = M.GetValueWrap(a, b, i-dx, j-dy, k-dz);

    // Centre point.
    double uc = M.GetValue(a, b, i, j, k);

    // Forward different approximation to first derivative.
    return (uc - ub) / h;
}

#endif // FINITEDIFFERENCEDEF