#ifndef ABSTRACTBDRYRULEBOUNCEBACKDEF
#define ABSTRACTBDRYRULEBOUNCEBACKDEF

#include "boundary/AbstractBoundaryRule.h"
#include "macroscopic/MacroscopicVariable.h"

template <class T>
class AbstractBdryRuleBounceBack : public AbstractBoundaryRule<T>
{
public:
    AbstractBdryRuleBounceBack() = default;

    AbstractBdryRuleBounceBack(AbstractLattice<T>* pDistribution)
        : AbstractBoundaryRule<T>(pDistribution) {}

    AbstractBdryRuleBounceBack(AbstractLattice<T>* pDistribution, MacroscopicVariable<T>* pDensity, T velx, T vely, T velz)
        : AbstractBoundaryRule<T>(pDistribution), mpDensity(pDensity), mWallVelX(velx), mWallVelY(vely), mWallVelZ(velz) {}

    ~AbstractBdryRuleBounceBack() = default;

    inline void SetWallVelocityX(T velocity)
    {
        mWallVelX = velocity;
    }

    inline void SetWallVelocityY(T velocity)
    {
        mWallVelY = velocity;
    }

    inline void SetWallVelocityZ(T velocity)
    {
        mWallVelZ = velocity;
    }

    void SetWallVelocityVector(T velx, T vely, T velz = 0)
    {
        mWallVelX = velz;
        mWallVelY = vely;
        mWallVelZ = velz;
    }

    void SetDensityPointer(MacroscopicVariable<T>* pDensity)
    {
        mpDensity = pDensity;
    }

protected:
    MacroscopicVariable<T>* mpDensity;

    T mWallVelX, mWallVelY, mWallVelZ;

    T compute_bounceback(int q, int i, int j, int k) const
    {
        // These lines are the same for all bounceback implementations.
        int qrev = (this->mpDistribution)->QRev(q);
        T u_dot_c = mWallVelX * (this->mpDistribution)->CX(qrev) + mWallVelY * (this->mpDistribution)->CY(qrev) + mWallVelZ * (this->mpDistribution)->CZ(qrev);
        T delta = static_cast<T>(6.0) * (this->mpDistribution)->W(q) * mpDensity->GetValue(i, j, k) * u_dot_c;
        return (this->mpDistribution)->GetPrevFStar(qrev, i, j, k) - delta;
    }
};

#endif // ABSTRACTBDRYRULEBOUNCEBACKDEF