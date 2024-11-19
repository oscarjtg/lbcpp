#ifndef FLUIDEVOLVERSRTDEF
#define FLUIDEVOLVERSRTDEF

#include "evolver/AbstractFluidEvolver.h"

template <typename T, int ND, int NQ>
class FluidEvolverSRT : public AbstractFluidEvolver<T, ND, NQ>
{
public:
    FluidEvolverSRT() = default;

    ~FluidEvolverSRT() = default;

    virtual void DoLocalCollision(AbstractLattice<T, ND, NQ>& f, std::array<T, NQ> flocal, T r_, T u_, T v_, T w_, T Fx_, T Fy_, T Fz_, int i, int j, int k) override;
};

#include "evolver/FluidEvolverSRT.tpp"

#endif // FLUIDEVOLVERSRTDEF