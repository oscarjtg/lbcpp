#ifndef LINEARFORCEANOMALYTEMP
#define LINEARFORCEANOMALYTEMP

#include "macroscopic/MacroscopicVariable.h"

template <class T>
class LinearForceAnomalyTemp
{
public:
    LinearForceAnomalyTemp(T ag_x,
                           T ag_y,
                           T ag_z,
                           T ref_temp,
                           int nx,
                           int ny,
                           int nz)
        : mAlphaGravityX(ag_x), mAlphaGravityY(ag_y), mAlphaGravityZ(ag_z),
        mReferenceTemperature(ref_temp), mNX(nx), mNY(ny), mNZ(nz) {}

    ~LinearForceAnomalyTemp() = default;

    // Expects ag_z < 0 for downwards gravity.
    void UpdateForce(MacroscopicVariable<T>& forceX,
                     MacroscopicVariable<T>& forceY,
                     MacroscopicVariable<T>& forceZ,
                     const MacroscopicVariable<T>& temperature)
    {
        T fx, fy, fz;
        for (int k = 0; k < mNZ; ++k)
        {
            for (int j = 0; j < mNY; ++j)
            {
                for (int i = 0; i < mNX; ++i)
                {
                    fx = -mAlphaGravityX * (temperature.GetValue(i, j, k) - mReferenceTemperature);
                    fy = -mAlphaGravityY * (temperature.GetValue(i, j, k) - mReferenceTemperature);
                    fz = -mAlphaGravityZ * (temperature.GetValue(i, j, k) - mReferenceTemperature);
                    forceX.SetValue(fx, i, j, k);
                    forceY.SetValue(fy, i, j, k);
                    forceZ.SetValue(fz, i, j, k);
                }
            }
        }
    }

    // Expects ag_z < 0 for downwards gravity.
    void UpdateForce(VectorField<T>& force,
                     const ScalarField<T>& temperature)
    {
        T fx, fy, fz;
        for (int k = 0; k < mNZ; ++k)
        {
            for (int j = 0; j < mNY; ++j)
            {
                for (int i = 0; i < mNX; ++i)
                {
                    fx = -mAlphaGravityX * (temperature.GetValue(i, j, k) - mReferenceTemperature);
                    fy = -mAlphaGravityY * (temperature.GetValue(i, j, k) - mReferenceTemperature);
                    fz = -mAlphaGravityZ * (temperature.GetValue(i, j, k) - mReferenceTemperature);
                    force.SetValue(fx, 0, i, j, k);
                    force.SetValue(fy, 1, i, j, k);
                    force.SetValue(fz, 2, i, j, k);
                }
            }
        }
    }

private:
    T mAlphaGravityX, mAlphaGravityY, mAlphaGravityZ, mReferenceTemperature;

    int mNX, mNY, mNZ;
};


#endif // LINEARFORCEANOMALYTEMP