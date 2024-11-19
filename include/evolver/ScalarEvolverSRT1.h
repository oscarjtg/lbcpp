#ifndef ScalarEvolverSRT1DEF
#define ScalarEvolverSRT1DEF

#include "evolver/AbstractScalarEvolver.h"

template <class T, int ND, int NQ>
class ScalarEvolverSRT1 : public AbstractScalarEvolver<T, ND, NQ>
{
public:
    ScalarEvolverSRT1() = default;

    ~ScalarEvolverSRT1() = default;
   
    virtual void DoLocalCollision(AbstractLattice<T, ND, NQ>& g, std::array<T, NQ> glocal, T c_, T u_, T v_, T w_, int i, int j, int k) override;
};

#include "evolver/ScalarEvolverSRT1.tpp"

#endif // ScalarEvolverSRT1DEF