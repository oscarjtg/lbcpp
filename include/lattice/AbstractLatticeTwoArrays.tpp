#include <iostream>
#include <string>
#include <vector>

#include "lattice/AbstractLatticeTwoArrays.h"

template <class T>
AbstractLatticeTwoArrays<T>::AbstractLatticeTwoArrays(int size_x, int size_y, int size_z, int dims, int vels, std::string id_string, std::string run_id, std::string savepath, int proc_num)
    : AbstractLattice<T>(size_x, size_y, size_z, size_x*size_y*size_z, dims, vels, id_string, run_id, savepath, proc_num), mFCurr(vels*size_x*size_y*size_z), mFPrev(vels*size_x*size_y*size_z) {}
    
template <class T>
void AbstractLatticeTwoArrays<T>::StreamDistributions()
{
    mFPrev.swap(mFCurr); // Swaps pointers internally.
}

template <class T>
void AbstractLatticeTwoArrays<T>::WriteToTextFile(const int timestep) const
{
    std::string basepath = this->construct_basepath(timestep);
    for (int q = 0; q < this->mVels; ++q)
    {
        this->write_txt(basepath, static_cast<T (AbstractLattice<T>::*)(int, int, int, int) const>(&AbstractLatticeTwoArrays<T>::GetCurrF), q);
    }
}
