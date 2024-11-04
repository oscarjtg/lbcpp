#ifndef NODEINFODEF
#define NODEINFODEF

#include <stdexcept>
#include <vector>

class NodeInfo
{
public:
    NodeInfo(int size_x, int size_y, int size_z)
        : mSizeX(size_x), mSizeY(size_y), mSizeZ(size_z), mCellTypeArray(size_x*size_y*size_z), mFillFraction(size_x*size_y*size_z) 
        {
            for (int id = 0; id < size_x*size_y*size_z; ++id)
            {
                mCellTypeArray[id] = mFluid;
                mFillFraction[id] = 1.0;
            }
        }

    ~NodeInfo() = default;

    inline void SetFluid(int i, int j, int k)
    {
        mCellTypeArray[idx(i, j, k)] = mFluid;
    }

    inline void SetSolid(int i, int j, int k)
    {
        mCellTypeArray[idx(i, j, k)] = mSolid;
    }

    inline void SetBoundary(int bdry_id, int i, int j, int k)
    {
        mCellTypeArray[idx(i, j, k)] = mBoundary + bdry_id;
    }

    void SetBoundaryOnBottom(int bdry_id)
    {
        int z_bot = 0;
        for (int j = 0; j < mSizeY; ++j)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                mCellTypeArray[idx(i, j, z_bot)] = mBoundary + bdry_id;
            }
        }
    }

    void SetBoundaryOnTop(int bdry_id)
    {
        int z_top = mSizeZ - 1;
        for (int j = 0; j < mSizeY; ++j)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                mCellTypeArray[idx(i, j, z_top)] = mBoundary + bdry_id;
            }
        }
    }

    void SetBoundaryOnLeft(int bdry_id)
    {
        int x_lef = 0;
        for (int k = 0; k < mSizeZ; ++k)
        {
            for (int j = 0; j < mSizeY; ++j)
            {
                mCellTypeArray[idx(x_lef, j, k)] = mBoundary + bdry_id;
            }
        }
    }

    void SetBoundaryOnRight(int bdry_id)
    {
        int x_rig = mSizeX - 1;
        for (int k = 0; k < mSizeZ; ++k)
        {
            for (int j = 0; j < mSizeY; ++j)
            {
                mCellTypeArray[idx(x_rig, j, k)] = mBoundary + bdry_id;
            }
        }
    }

    void SetBoundaryOnFront(int bdry_id)
    {
        int y_fro = 0;
        for (int k = 0; k < mSizeZ; ++k)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                mCellTypeArray[idx(i, y_fro, k)] = mBoundary + bdry_id;
            }
        }
    }

    void SetBoundaryOnBack(int bdry_id)
    {
        int y_bac = mSizeY - 1;
        for (int k = 0; k < mSizeZ; ++k)
        {
            for (int i = 0; i < mSizeX; ++i)
            {
                mCellTypeArray[idx(i, y_bac, k)] = mBoundary + bdry_id;
            }
        }
    }

    int GetBoundaryID(int i, int j, int k) const
    {
        int bdry_id = mCellTypeArray[idx(i, j, k)] - mBoundary;
        if (bdry_id < 0)
        {
            throw std::runtime_error("Boundary ID is negative.");
        }
        else
        {
            return bdry_id;
        }
    }

    inline void SetInterface(int i, int j, int k)
    {
        mCellTypeArray[idx(i, j, k)] = mInterface;
    }

    inline void SetGas(int i, int j, int k)
    {
        mCellTypeArray[idx(i, j, k)] = mGas;
    }

    inline bool IsFluid(int i, int j, int k) const
    {
        return mCellTypeArray[idx(i, j, k)] == mFluid;
    }

    inline bool IsSolid(int i, int j, int k) const
    {
        return mCellTypeArray[idx(i, j, k)] == mSolid;
    }

    inline bool IsBoundary(int i, int j, int k) const
    {
        return mCellTypeArray[idx(i, j, k)] >= mBoundary;
    }

    inline bool IsInterface(int i, int j, int k) const
    {
        return mCellTypeArray[idx(i, j, k)] == mInterface;
    }

    inline bool IsGas(int i, int j, int k) const
    {
        return mCellTypeArray[idx(i, j, k)] == mGas;
    }

    inline int GetNX() const
    {
        return mSizeX;
    }

    inline int GetNY() const
    {
        return mSizeY;
    }

    inline int GetNZ() const
    {
        return mSizeZ;
    }

private:
    int mSizeX, mSizeY, mSizeZ;

    const int mSolid = 0, mFluid = 1, mGas = 2, mInterface = 3, mBoundary = 4;

    const double tol = 1.0e-3;

    std::vector<int> mCellTypeArray;

    std::vector<double> mFillFraction;

    inline int idx(int i, int j, int k) const
    {
        return i + mSizeX * j + mSizeX * mSizeY * k;
    }
};

#endif // NODEINFODEF