#include <iostream>
#include <string>
#include <vector>

#include "lattice/AbstractLatticeTwoArrays.h"

template <typename T, int ND, int NQ>
AbstractLatticeTwoArrays<T, ND, NQ>::AbstractLatticeTwoArrays(int size_x, int size_y, int size_z, std::string id_string, std::string run_id, std::string savepath, int proc_num)
    : AbstractLattice<T, ND, NQ>(size_x, size_y, size_z, size_x*size_y*size_z, id_string, run_id, savepath, proc_num), mFCurr(NQ*size_x*size_y*size_z), mFPrev(NQ*size_x*size_y*size_z) {}
    
template <typename T, int ND, int NQ>
void AbstractLatticeTwoArrays<T, ND, NQ>::StreamDistributions()
{
    mFPrev.swap(mFCurr); // Swaps pointers internally.
}

template <typename T, int ND, int NQ>
void AbstractLatticeTwoArrays<T, ND, NQ>::WriteToTextFile(const int timestep) const
{
    std::string basepath = this->construct_basepath(timestep);
    for (int q = 0; q < NQ; ++q)
    {
        this->write_txt(basepath, static_cast<T (AbstractLattice<T, ND, NQ>::*)(int, int, int, int) const>(&AbstractLatticeTwoArrays<T, ND, NQ>::GetCurrF), q);
    }
}
