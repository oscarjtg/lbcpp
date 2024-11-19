/*
DO NOT USE -- BLOWS UP!
*/

#ifndef FluidEvolverSRT3DEF
#define FluidEvolverSRT3DEF

#include "evolver/AbstractFluidEvolver.h"

template <class T, int ND, int NQ>
class FluidEvolverSRT3 : public AbstractFluidEvolver<T, ND, NQ>
{
public:
    FluidEvolverSRT3() = default;

    ~FluidEvolverSRT3() = default;

    virtual void DoLocalCollision(AbstractLattice<T, ND, NQ>& f, std::array<T, NQ> flocal, T r_, T u_, T v_, T w_, T Fx_, T Fy_, T Fz_, int i, int j, int k) override;
};

#include "evolver/FluidEvolverSRT3.tpp"

#endif // FluidEvolverSRT3DEF