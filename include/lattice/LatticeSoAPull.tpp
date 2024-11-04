#include <iostream>
#include <string>

#include "lattice/LatticeSoAPull.h"

template <class T>
LatticeSoAPull<T>::LatticeSoAPull(int size_x, int size_y, int size_z, int dims, int vels, std::string id_string, std::string run_id, std::string savepath, int proc_num)
    : AbstractLatticeTwoArrays<T>(size_x, size_y, size_z, dims, vels, id_string, run_id, savepath, proc_num) {}
    
template <class T>
T LatticeSoAPull<T>::GetPrevFStar(int q, int i, int j, int k) const
{
    return this->mFPrev[idx(q, i, j, k)];
}

template <class T>
T LatticeSoAPull<T>::GetCurrF(int q, int i, int j, int k) const
{
    int ii = (i - this->mEX[q] + this->mSizeX) % this->mSizeX;
    int jj = (j - this->mEY[q] + this->mSizeY) % this->mSizeY;
    int kk = (k - this->mEZ[q] + this->mSizeZ) % this->mSizeZ;

    return this->mFPrev[idx(q, ii, jj, kk)];
}

template <class T>
T LatticeSoAPull<T>::GetCurrFStar(int q, int i, int j, int k) const
{
    return this->mFCurr[idx(q, i, j, k)];
}

template <class T>
void LatticeSoAPull<T>::SetCurrF(T f_, int q, int i, int j, int k)
{
    int ii = (i - this->mEX[q] + this->mSizeX) % this->mSizeX;
    int jj = (j - this->mEY[q] + this->mSizeY) % this->mSizeY;
    int kk = (k - this->mEZ[q] + this->mSizeZ) % this->mSizeZ;

    this->mFPrev[idx(q, ii, jj, kk)] = f_;
}

template <class T>
void LatticeSoAPull<T>::SetCurrFStar(T f_, int q, int i, int j, int k)
{
    this->mFCurr[idx(q, i, j, k)] = f_;
}

template <class T>
void LatticeSoAPull<T>::WriteToTextFile(const int timestep) const
{   
    std::string basepath = this->construct_basepath(timestep);
    for (int q = 0; q < this->mVels; ++q)
    {
        this->write_txt(basepath, static_cast<T (AbstractLattice<T>::*)(int, int, int, int) const>(&LatticeSoAPull<T>::GetCurrF), q);
    }
}

