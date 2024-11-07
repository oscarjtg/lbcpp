#ifndef DIAGNOSTICNUSSELTDEF
#define DIAGNOSTICNUSSELTDEF

#include <iostream>

#include "macroscopic/MacroscopicVariable.h"

template <class T>
class DiagnosticNusselt
{
public:
    DiagnosticNusselt(int nx, int ny, int nz, double kappa) : mNX(nx), mNY(ny), mNZ(nz), mDiffusivity(kappa) {}

    ~DiagnosticNusselt() = default;

    double ComputeNusseltNumber(const MacroscopicVariable<T>& velz, const MacroscopicVariable<T>& temp)
    {
        // First compute average of velz * temp over whole domain.
        double wT = 0.0;
        for (int k = 0; k < mNZ; ++k)
        {
            for (int j = 0; j < mNY; ++j)
            {
                for (int i = 0; i < mNX; ++i)
                {
                    wT += static_cast<double>(velz.GetValue(i, j, k)) * static_cast<double>(temp.GetValue(i, j, k));
                }
            }
        }
        wT /= (mNX * mNY); // This is the averaging part. 
        // Note did not divide by mNZ because would have to multiply by it again in next line for Nu.
        double nusselt_number = 1.0 + wT / mDiffusivity;
        std::cout << "Nusselt number = " << nusselt_number << "\n";
        return nusselt_number;
    }

private:
    int mNX, mNY, mNZ;
    double mDiffusivity;
};

#endif // DIAGNOSTICNUSSELTDEF