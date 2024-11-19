#ifndef ABSTRACTLATTICETWOARRAYSDEF
#define ABSTRACTLATTICETWOARRAYSDEF

#include "lattice/AbstractLattice.h"
#include <vector>
#include <string>

template <typename T, int ND, int NQ>
class AbstractLatticeTwoArrays : public AbstractLattice<T, ND, NQ>
{
public:
    AbstractLatticeTwoArrays(int size_x, int size_y, int size_z, std::string id_string = "f", std::string run_id = "test", std::string savepath = "output", int proc_num = 0);

    ~AbstractLatticeTwoArrays() = default;

    void StreamDistributions() override;

    void WriteToTextFile(const int timestep = 0) const override;

protected:
    std::vector<T> mFCurr;
    std::vector<T> mFPrev;
};

#include "lattice/AbstractLatticeTwoArrays.tpp"

#endif // ABSTRACTLATTICETWOARRAYSDEF