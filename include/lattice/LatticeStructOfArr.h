#ifndef LATTICESTRUCTOFARRDEF
#define LATTICESTRUCTOFARRDEF

#include "lattice/AbstractLattice.h"
#include <vector>
#include <string>

template <class T>
class LatticeStructOfArr : public AbstractLattice<T> 
{
public:
    LatticeStructOfArr(int size_x, int size_y, int size_z, int dims, int vels, std::string id_string = "f", std::string run_id = "test", std::string savepath = "output", int proc_num = 0);

    ~LatticeStructOfArr() = default;

    T GetPrevFStar(int q, int i, int j, int k) const override;

    T GetCurrF(int q, int i, int j, int k) const override;

    T GetCurrFStar(int q, int i, int j, int k) const override;

    void SetCurrF(T f_, int q, int i, int j, int k) override;

    void SetCurrFStar(T f_, int q, int i, int j, int k) override;

    void StreamDistributions() override;

    void WriteToTextFile(const int timestep = 0) const override;

private:
    std::vector<T> mFCurr;
    std::vector<T> mFPrev;

    inline int idx(int q, int i, int j, int k) const
    {
        return i + this->mSizeX * j + this->mSizeX * this->mSizeY * k + this->mSizeX * this->mSizeY * this->mSizeZ * q;
    }
};

#endif // LATTICESTRUCTOFARRDEF