#ifndef ABSTRACTLATTICEDEF
#define ABSTRACTLATTICEDEF

#include <array>
#include <string>

template <class T>
class AbstractLattice 
{
public:
    AbstractLattice(int size_x, int size_y, int size_z, int grid_size, int dims, int vels, std::string id_string = "f", std::string run_id = "test", std::string savepath = "output", int proc_num = 0);

    virtual ~AbstractLattice() = default;

    virtual void DisplayLatticeParameters() const;

    virtual T GetPrevFStar(int q, int i, int j, int k) const = 0;

    virtual T GetCurrF(int q, int i, int j, int k) const = 0;

    virtual T GetCurrFStar(int q, int i, int j, int k) const = 0;

    virtual void SetCurrF(T f_, int q, int i, int j, int k) = 0;

    virtual void SetCurrFStar(T f_, int q, int i, int j, int k) = 0;

    virtual void StreamDistributions() = 0;

    virtual void WriteToTextFile(const int timestep = 0) const = 0;

    inline int EX(int q) const
    {
        return mEX[q];
    }

    inline int EY(int q) const
    {
        return mEY[q];
    }

    inline int EZ(int q) const
    {
        return mEZ[q];
    }

    inline T CX(int q) const
    {
        return static_cast<T>(mEX[q]);
    }

    inline T CY(int q) const
    {
        return static_cast<T>(mEY[q]);
    }

    inline T CZ(int q) const
    {
        return static_cast<T>(mEZ[q]);
    }    

    inline T W(int q) const
    {
        return mW[q];
    }

    inline int QRev(int q) const
    {
        return mQRev[q];
    }

    inline int ND() const
    {
        return mDims;
    }

    inline int NQ() const
    {
        return mVels;
    }

    inline T CSI() const
    {
        return mCSI;
    }

protected:
    const int mSizeX, mSizeY, mSizeZ, mGridSize, mDims, mVels;
    
    const std::string mIDString, mRunID, mSaveDirectory;
    const int mProcessNumber;
    
    std::array<int, 27> mEX;
    std::array<int, 27> mEY;
    std::array<int, 27> mEZ;
    std::array<T, 27> mW;
    std::array<int, 27> mQRev;

    T mCSI; // Lattice speed of sound, inverse, squared.

    std::string construct_basepath(const int timestep) const;

    void write_txt(const std::string& path, T (AbstractLattice<T>::*get_func)(int, int, int, int) const, int q) const;
};

#include "lattice/AbstractLattice.tpp"

#endif // ABSTRACTLATTICEDEF