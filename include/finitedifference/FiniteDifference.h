#include "macroscopic/MacroscopicVariable.h"

#include <cmath>

template <typename T>
double ComputeCentredDiff(const MacroscopicVariable<T>& M, int i, int j, int k, int dx, int dy, int dz)
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
double ComputeForwardDiff(const MacroscopicVariable<T>& M, int i, int j, int k, int dx, int dy, int dz)
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
double ComputeBackwardDiff(const MacroscopicVariable<T>& M, int i, int j, int k, int dx, int dy, int dz)
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