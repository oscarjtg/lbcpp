#ifndef SCALAREVOLVERSRTDEF
#define SCALAREVOLVERSRTDEF

#include "evolver/AbstractScalarEvolver.h"

template <class T>
class ScalarEvolverSRT : public AbstractScalarEvolver<T>
{
public:
    ScalarEvolverSRT() = default;

    ~ScalarEvolverSRT() = default;

    void SetScalarDiffusivity(T kappa) override;
   
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

#include "evolver/ScalarEvolverSRT.tpp"

#endif // SCALAREVOLVERSRTDEF