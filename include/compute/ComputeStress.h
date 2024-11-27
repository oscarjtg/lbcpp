#ifndef COMPUTESTRESSDEF
#define COMPUTESTRESSDEF

#include "macroscopic/AllFields.h"
#include "finitedifference/FiniteDifference.h"

template <typename T>
void ComputeStressFiniteDifferencePeriodic(TensorField<T>& sigma, const VectorField<T>& vel, const double nu, const double rho0)
{
    static constexpr std::array<int, 3> DX = {1, 0, 0};
    static constexpr std::array<int, 3> DY = {0, 1, 0};
    static constexpr std::array<int, 3> DZ = {0, 0, 1};

    for (int k = 0; k < sigma.GetNZ(); ++k)
    {
        for (int j = 0; j < sigma.GetNY(); ++j)
        {
            for (int i = 0; i < sigma.GetNZ(); ++i)
            {
                for (int a = 0; a < 3; ++a)
                {
                    for (int b = 0; b < 3; ++b)
                    {
                        double e1 = ComputeCentredDiff(vel, a, i, j, k, DX[b], DY[b], DZ[b]);
                        double e2 = ComputeCentredDiff(vel, b, i, j, k, DX[a], DY[a], DZ[a]);
                        double e_ = e1 + e2;
                        double s_ = rho0 * nu * e_;
                        sigma.SetValue(s_, a, b, i, j, k);
                    }
                }
            }
        }
    }
}

#endif // COMPUTESTRESSDEF