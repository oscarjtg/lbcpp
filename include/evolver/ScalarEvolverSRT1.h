#ifndef ScalarEvolverSRT1DEF
#define ScalarEvolverSRT1DEF

#include "evolver/AbstractScalarEvolver.h"

template <class T>
class ScalarEvolverSRT1 : public AbstractScalarEvolver<T>
{
public:
    ScalarEvolverSRT1() = default;

    ~ScalarEvolverSRT1() = default;

    void SetScalarDiffusivity(AbstractLattice<T>& g, T kappa) override;
   
    void Initialise(AbstractLattice<T>& g,
                    const MacroscopicVariable<T>& conc,
                    const MacroscopicVariable<T>& velx,
                    const MacroscopicVariable<T>& vely,
                    const MacroscopicVariable<T>& velz,
                    const NodeInfo& node,
                    const BoundaryInfo<T>& bdry) override;
    
    void DoTimestep(AbstractLattice<T>& g,
                    MacroscopicVariable<T>& conc,
                    const MacroscopicVariable<T>& velx,
                    const MacroscopicVariable<T>& vely,
                    const MacroscopicVariable<T>& velz,
                    const NodeInfo& node,
                    const BoundaryInfo<T>& bdry) override;

private:
    T mOmega;
};

#include "evolver/ScalarEvolverSRT1.tpp"

#endif // ScalarEvolverSRT1DEF