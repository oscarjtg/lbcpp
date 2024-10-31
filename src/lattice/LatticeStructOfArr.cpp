#include <iostream>
#include <string>

#include "lattice/LatticeStructOfArr.h"

template <class T>
LatticeStructOfArr<T>::LatticeStructOfArr(int size_x, int size_y, int size_z, int dims, int vels, std::string id_string, std::string run_id, std::string savepath, int proc_num)
    : AbstractLattice<T>(size_x, size_y, size_z, size_x*size_y*size_z, dims, vels, id_string, run_id, savepath, proc_num), mFCurr(size_x*size_y*size_z), mFPrev(size_x*size_y*size_z) {}
    
template <class T>
T LatticeStructOfArr<T>::GetPrevFStar(int q, int i, int j, int k) const
{
    return mFPrev[idx(q, i, j, k)];
}

template <class T>
T LatticeStructOfArr<T>::GetCurrF(int q, int i, int j, int k) const
{
    int ii = (i - this->mEX[q] + this->mSizeX) % this->mSizeX;
    int jj = (j - this->mEY[q] + this->mSizeY) % this->mSizeY;
    int kk = (k - this->mEZ[q] + this->mSizeZ) % this->mSizeZ;

    return mFPrev[idx(q, ii, jj, kk)];
}

template <class T>
T LatticeStructOfArr<T>::GetCurrFStar(int q, int i, int j, int k) const
{
    return mFCurr[idx(q, i, j, k)];
}

template <class T>
void LatticeStructOfArr<T>::SetCurrF(T f_, int q, int i, int j, int k)
{
    int ii = (i - this->mEX[q] + this->mSizeX) % this->mSizeX;
    int jj = (j - this->mEY[q] + this->mSizeY) % this->mSizeY;
    int kk = (k - this->mEZ[q] + this->mSizeZ) % this->mSizeZ;

    mFPrev[idx(q, ii, jj, kk)] = f_;
}

template <class T>
void LatticeStructOfArr<T>::SetCurrFStar(T f_, int q, int i, int j, int k)
{
    mFCurr[idx(q, i, j, k)] = f_;
}

template <class T>
void LatticeStructOfArr<T>::StreamDistributions()
{
    mFPrev.swap(mFCurr); // Swaps pointers internally.
}

template <class T>
void LatticeStructOfArr<T>::WriteToTextFile(const int timestep) const
{
    std::string basepath = this->construct_basepath(timestep);
    for (int q = 0; q < this->mVels; ++q)
    {
        this->write_txt(basepath, static_cast<T (AbstractLattice<T>::*)(int, int, int, int) const>(&LatticeStructOfArr<T>::GetCurrF), q);
    }
}

// Explicit template instantiation.
//template class LatticeStructOfArr<_Float16>;
template class LatticeStructOfArr<float>;
template class LatticeStructOfArr<double>;