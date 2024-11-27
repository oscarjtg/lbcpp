#ifndef ABSTRACTBDRYRULEANTIBOUNCEBACKDEF
#define ABSTRACTBDRYRULEANTIBOUNCEBACKDEF

#include "boundary/AbstractBoundaryRule.h"
#include "equilibria/Equilibria.h"
#include "macroscopic/AllFields.h"

template <typename T, int ND, int NQ>
class AbstractBdryRuleAntiBounceBack : public AbstractBoundaryRule<T, ND, NQ>
{
public:
    AbstractBdryRuleAntiBounceBack() = default;

    AbstractBdryRuleAntiBounceBack(AbstractLattice<T, ND, NQ>* pDistribution)
        : AbstractBoundaryRule<T, ND, NQ>(pDistribution) {}

    AbstractBdryRuleAntiBounceBack(AbstractLattice<T, ND, NQ>* pDistribution, T wall_conc, T velx, T vely, T velz)
        : AbstractBoundaryRule<T, ND, NQ>(pDistribution), mWallConc(wall_conc), mWallVelX(velx), mWallVelY(vely), mWallVelZ(velz) {}

    ~AbstractBdryRuleAntiBounceBack() = default;

    inline void SetWallConcentration(T wall_conc)
    {
        mWallConc = wall_conc;
    }

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

protected:
    T mWallConc, mWallVelX, mWallVelY, mWallVelZ;

    T compute_antibounceback(int q, int i, int j, int k) const
    {
        // These lines are the same for all bounceback implementations.
        int qrev = (this->mpDistribution)->QRev(q);
        T weight = (this->mpDistribution)->W(q);
        T u_dot_c = mWallVelX * (this->mpDistribution)->CX(qrev) + mWallVelY * (this->mpDistribution)->CY(qrev) + mWallVelZ * (this->mpDistribution)->CZ(qrev);
        T usq = mWallVelX*mWallVelX + mWallVelY*mWallVelY + mWallVelZ*mWallVelZ;
        //std::cout << "compute_antibounceback(q, i, j, k): ";
        //std::cout << "qrev = " << qrev;
        //std::cout << ", mWallConc = " << mWallConc;
        //std::cout << ", weight = " << weight;
        //std::cout << ", u_dot_c = " << u_dot_c;
        //std::cout << ", usq = " << usq << "\n";
        return -(this->mpDistribution)->GetPrevFStar(qrev, i, j, k) + static_cast<T>(2.0)*computeSecondOrderEquilibrium(mWallConc, u_dot_c, usq, weight);
    }
};

#endif // ABSTRACTBDRYRULEANTIBOUNCEBACKDEF