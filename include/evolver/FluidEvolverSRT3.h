#ifndef FluidEvolverSRT3DEF
#define FluidEvolverSRT3DEF

#include "evolver/AbstractFluidEvolver.h"

template <class T>
class FluidEvolverSRT3 : public AbstractFluidEvolver<T>
{
public:
    FluidEvolverSRT3() = default;

    ~FluidEvolverSRT3() = default;

    void SetKinematicViscosity(AbstractLattice<T>& f, T nu) override;

    void SetBulkViscosity(AbstractLattice<T>& f, T eta) override;

    void Initialise(AbstractLattice<T>& f,
                    const MacroscopicVariable<T>& dens,
                    const MacroscopicVariable<T>& velx,
                    const MacroscopicVariable<T>& vely,
                    const MacroscopicVariable<T>& velz,
                    const MacroscopicVariable<T>& Fx,
                    const MacroscopicVariable<T>& Fy,
                    const MacroscopicVariable<T>& Fz,
                    const NodeInfo& node,
                    const BoundaryInfo<T>& bdry) override;

    void DoTimestep(AbstractLattice<T>& f,
                    MacroscopicVariable<T>& dens,
                    MacroscopicVariable<T>& velx,
                    MacroscopicVariable<T>& vely,
                    MacroscopicVariable<T>& velz,
                    const MacroscopicVariable<T>& Fx,
                    const MacroscopicVariable<T>& Fy,
                    const MacroscopicVariable<T>& Fz,
                    NodeInfo& node,
                    BoundaryInfo<T>& bdry) override;

private:
    T mOmega;
};

#include "evolver/FluidEvolverSRT3.tpp"

#endif // FluidEvolverSRT3DEF