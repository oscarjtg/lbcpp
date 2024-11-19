#include <iostream>
#include <string>

#include "lattice/LatticeSoAPull.h"

template <typename T, int ND, int NQ>
LatticeSoAPull<T, ND, NQ>::LatticeSoAPull(int size_x, int size_y, int size_z, std::string id_string, std::string run_id, std::string savepath, int proc_num)
    : AbstractLatticeTwoArrays<T, ND, NQ>(size_x, size_y, size_z, id_string, run_id, savepath, proc_num) {}
    
template <typename T, int ND, int NQ>
T LatticeSoAPull<T, ND, NQ>::GetPrevFStar(int q, int i, int j, int k) const
{
    return this->mFPrev[idx(q, i, j, k)];
}

template <typename T, int ND, int NQ>
T LatticeSoAPull<T, ND, NQ>::GetCurrF(int q, int i, int j, int k) const
{
    int ii = (i - this->mLP.mEX[q] + this->mSizeX) % this->mSizeX;
    int jj = (j - this->mLP.mEY[q] + this->mSizeY) % this->mSizeY;
    int kk = (k - this->mLP.mEZ[q] + this->mSizeZ) % this->mSizeZ;

    return this->mFPrev[idx(q, ii, jj, kk)];
}

template <typename T, int ND, int NQ>
T LatticeSoAPull<T, ND, NQ>::GetCurrFStar(int q, int i, int j, int k) const
{
    return this->mFCurr[idx(q, i, j, k)];
}

template <typename T, int ND, int NQ>
void LatticeSoAPull<T, ND, NQ>::SetCurrF(T f_, int q, int i, int j, int k)
{
    int ii = (i - this->mLP.mEX[q] + this->mSizeX) % this->mSizeX;
    int jj = (j - this->mLP.mEY[q] + this->mSizeY) % this->mSizeY;
    int kk = (k - this->mLP.mEZ[q] + this->mSizeZ) % this->mSizeZ;

    this->mFPrev[idx(q, ii, jj, kk)] = f_;
}

template <typename T, int ND, int NQ>
void LatticeSoAPull<T, ND, NQ>::SetCurrFStar(T f_, int q, int i, int j, int k)
{
    this->mFCurr[idx(q, i, j, k)] = f_;
}

template <typename T, int ND, int NQ>
void LatticeSoAPull<T, ND, NQ>::WriteToTextFile(const int timestep) const
{   
    std::string basepath = this->construct_basepath(timestep);
    for (int q = 0; q < NQ; ++q)
    {
        this->write_txt(basepath, static_cast<T (AbstractLattice<T, ND, NQ>::*)(int, int, int, int) const>(&LatticeSoAPull<T, ND, NQ>::GetCurrF), q);
    }
}

