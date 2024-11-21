#ifndef ABSTRACTLATTICEDEF
#define ABSTRACTLATTICEDEF

#include <array>
#include <string>

template <typename T, int ND, int NQ>
struct LatticeParameters {};

// Specialisation for D2Q4
template <typename T>
struct LatticeParameters<T, 2, 4> : LatticeParameters<T, 0, 0>
{
    static constexpr std::array<int, 4> mEX = {1, 0, -1, 0};
    static constexpr std::array<int, 4> mEY = {0, 0, 0, 0};
    static constexpr std::array<int, 4> mEZ = {0, 1, 0, -1};
    static constexpr std::array<T, 4> mW = {1./4., 1./4., 1./4., 1./4.};
    static constexpr std::array<int, 4> mQRev = {2, 3, 0, 1};
    const T mCSI = 2;
};

// Specialisation for D2Q5
template <typename T>
struct LatticeParameters<T, 2, 5> : LatticeParameters<T, 0, 0>
{
    static constexpr std::array<int, 5> mEX = {0, 1, 0, -1, 0};
    static constexpr std::array<int, 5> mEY = {0, 0, 0, 0, 0};
    static constexpr std::array<int, 5> mEZ = {0, 0, 1, 0, -1};
    static constexpr std::array<T, 5> mW = {1./3., 1./6., 1./6., 1./6., 1./6.};
    static constexpr std::array<int, 5> mQRev = {0, 3, 4, 1, 2};
    const T mCSI = 3;
};

// Specialisation for D2Q9
template <typename T>
struct LatticeParameters<T, 2, 9> : LatticeParameters<T, 0, 0>
{
    static constexpr std::array<int, 9> mEX = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    static constexpr std::array<int, 9> mEY = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    static constexpr std::array<int, 9> mEZ = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    static constexpr std::array<T, 9> mW = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
    static constexpr std::array<int, 9> mQRev = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    const T mCSI = 3;
};

// Specialisation for D3Q6
template <typename T>
struct LatticeParameters<T, 3, 6> : LatticeParameters<T, 0, 0>
{
    static constexpr std::array<int, 6> mEX = {1, -1, 0, 0, 0, 0};
    static constexpr std::array<int, 6> mEY = {0, 0, 1, -1, 0, 0};
    static constexpr std::array<int, 6> mEZ = {0, 0, 0, 0, 1, -1};
    static constexpr std::array<T, 6> mW = {1./6., 1./6., 1./6., 1./6., 1./6., 1./6.};
    static constexpr std::array<int, 6> mQRev = {1, 0, 3, 2, 5, 4};
    const T mCSI = 2;
};

// Specialisation for D3Q7
template <typename T>
struct LatticeParameters<T, 3, 7> : LatticeParameters<T, 0, 0>
{
    static constexpr std::array<int, 7> mEX = {0, 1, -1, 0, 0, 0, 0};
    static constexpr std::array<int, 7> mEY = {0, 0, 0, 1, -1, 0, 0};
    static constexpr std::array<int, 7> mEZ = {0, 0, 0, 0, 0, 1, -1};
    static constexpr std::array<T, 7> mW = {1./4., 1./8., 1./8., 1./8., 1./8., 1./8., 1./8.};
    static constexpr std::array<int, 7> mQRev = {0, 2, 1, 4, 3, 6, 5};
    const T mCSI = 2;
};

// Specialisation for D3Q15
template <typename T>
struct LatticeParameters<T, 3, 15> : LatticeParameters<T, 0, 0>
{
    static constexpr std::array<int, 15> mEX = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1};
    static constexpr std::array<int, 15> mEY = {0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1};
    static constexpr std::array<int, 15> mEZ = {0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1};
    static constexpr std::array<T, 15> mW = {2./9., 1./9., 1./9., 1./9., 1./9., 1./9., 1./9., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72.};
    static constexpr std::array<int, 15> mQRev = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};
    const T mCSI = 3;
};

// Specialisation for D3Q19
template <typename T>
struct LatticeParameters<T, 3, 19> : LatticeParameters<T, 0, 0>
{
    static constexpr std::array<int, 19> mEX = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0};
    static constexpr std::array<int, 19> mEY = {0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1};
    static constexpr std::array<int, 19> mEZ = {0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1};
    static constexpr std::array<T, 19> mW = {1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,  1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.};
    static constexpr std::array<int, 19> mQRev = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};
    const T mCSI = 3;
};

// Specialisation for D3Q27
template <typename T>
struct LatticeParameters<T, 3, 27> : LatticeParameters<T, 0, 0>
{
    static constexpr std::array<int, 27> mEX = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1};
    static constexpr std::array<int, 27> mEY = {0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1};
    static constexpr std::array<int, 27> mEZ = {0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1};
    static constexpr std::array<T, 27> mW = {8./27., 2./27., 2./27., 2./27., 2./27., 2./27., 2./27., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216.};
    static constexpr std::array<int, 27> mQRev = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};
    const T mCSI = 3;
};


template <typename T, int ND, int NQ>
class AbstractLattice 
{
public:
    AbstractLattice(int size_x, int size_y, int size_z, int grid_size, std::string id_string = "f", std::string run_id = "test", std::string savepath = "output", int proc_num = 0);

    virtual ~AbstractLattice() = default;

    virtual void DisplayLatticeParameters() const;

    virtual T GetPrevFStar(int q, int i, int j, int k) const = 0;

    virtual T GetCurrF(int q, int i, int j, int k) const = 0;

    virtual T GetCurrFStar(int q, int i, int j, int k) const = 0;

    virtual void SetCurrF(T f_, int q, int i, int j, int k) = 0;

    virtual void SetCurrFStar(T f_, int q, int i, int j, int k) = 0;

    virtual void StreamDistributions() = 0;

    virtual void WriteToTextFile(const int timestep = 0) const = 0;

    // Getters.
    inline int EX(int q) const { return mLP.mEX[q]; }
    inline int EY(int q) const { return mLP.mEY[q]; }
    inline int EZ(int q) const { return mLP.mEZ[q]; }
    inline T CX(int q) const { return static_cast<T>(mLP.mEX[q]); }
    inline T CY(int q) const { return static_cast<T>(mLP.mEY[q]); }
    inline T CZ(int q) const { return static_cast<T>(mLP.mEZ[q]); }    
    inline T W(int q) const { return mLP.mW[q]; }
    inline int QRev(int q) const { return mLP.mQRev[q]; }
    inline T CSI() const { return mLP.mCSI; }
    inline int GetNX() const { return mSizeX; }
    inline int GetNY() const { return mSizeY; }
    inline int GetNZ() const { return mSizeZ; }

protected:
    const int mSizeX, mSizeY, mSizeZ, mGridSize;
    
    const std::string mIDString, mRunID, mSaveDirectory;
    const int mProcessNumber;
    
    const LatticeParameters<T, ND, NQ> mLP;

    std::string construct_basepath(const int timestep) const;

    void write_txt(const std::string& path, T (AbstractLattice<T, ND, NQ>::*get_func)(int, int, int, int) const, int q) const;
};

/*
// Specialisation for D2Q4
template <typename T>
class AbstractLattice<T, 2, 4> : public AbstractLattice<T, 0, 0>
{
public:
    inline int EX(int q) const { return mEX[q]; }
    inline int EY(int q) const { return mEY[q]; }
    inline int EZ(int q) const { return mEZ[q]; }
    inline T CX(int q) const { return static_cast<T>(mEX[q]); }
    inline T CY(int q) const { return static_cast<T>(mEY[q]); }
    inline T CZ(int q) const { return static_cast<T>(mEZ[q]); }    
    inline T W(int q) const { return mW[q]; }
    inline int QRev(int q) const { return mQRev[q]; }
    inline T CSI() const { return mCSI; }
protected:
    static constexpr std::array<int, 4> mEX = {1, 0, -1, 0};
    static constexpr std::array<int, 4> mEY = {0, 0, 0, 0};
    static constexpr std::array<int, 4> mEZ = {0, 1, 0, -1};
    static constexpr std::array<T, 4> mW = {1./4., 1./4., 1./4., 1./4.};
    static constexpr std::array<int, 4> mQRev = {2, 3, 0, 1};
    const T mCSI = 2;
};

// Specialisation for D2Q5
template <typename T>
class AbstractLattice<T, 2, 5> : public AbstractLattice<T, 0, 0>
{
public:
    inline int EX(int q) const { return mEX[q]; }
    inline int EY(int q) const { return mEY[q]; }
    inline int EZ(int q) const { return mEZ[q]; }
    inline T CX(int q) const { return static_cast<T>(mEX[q]); }
    inline T CY(int q) const { return static_cast<T>(mEY[q]); }
    inline T CZ(int q) const { return static_cast<T>(mEZ[q]); }    
    inline T W(int q) const { return mW[q]; }
    inline int QRev(int q) const { return mQRev[q]; }
    inline T CSI() const { return mCSI; }
protected:
    static constexpr std::array<int, 5> mEX = {0, 1, 0, -1, 0};
    static constexpr std::array<int, 5> mEY = {0, 0, 0, 0, 0};
    static constexpr std::array<int, 5> mEZ = {0, 0, 1, 0, -1};
    static constexpr std::array<T, 5> mW = {1./3., 1./6., 1./6., 1./6., 1./6.};
    static constexpr std::array<int, 5> mQRev = {0, 3, 4, 1, 2};
    const T mCSI = 3;
};

// Specialisation for D2Q9
template <typename T>
class AbstractLattice<T, 2, 9> : public AbstractLattice<T, 0, 0>
{
public:
    inline int EX(int q) const { return mEX[q]; }
    inline int EY(int q) const { return mEY[q]; }
    inline int EZ(int q) const { return mEZ[q]; }
    inline T CX(int q) const { return static_cast<T>(mEX[q]); }
    inline T CY(int q) const { return static_cast<T>(mEY[q]); }
    inline T CZ(int q) const { return static_cast<T>(mEZ[q]); }    
    inline T W(int q) const { return mW[q]; }
    inline int QRev(int q) const { return mQRev[q]; }
    inline T CSI() const { return mCSI; }
protected:
    static constexpr std::array<int, 9> mEX = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    static constexpr std::array<int, 9> mEY = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    static constexpr std::array<int, 9> mEZ = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    static constexpr std::array<T, 9> mW = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
    static constexpr std::array<int, 9> mQRev = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    const T mCSI = 3;
};

// Specialisation for D3Q6
template <typename T>
class AbstractLattice<T, 3, 6> : public AbstractLattice<T, 0, 0>
{
public:
    inline int EX(int q) const { return mEX[q]; }
    inline int EY(int q) const { return mEY[q]; }
    inline int EZ(int q) const { return mEZ[q]; }
    inline T CX(int q) const { return static_cast<T>(mEX[q]); }
    inline T CY(int q) const { return static_cast<T>(mEY[q]); }
    inline T CZ(int q) const { return static_cast<T>(mEZ[q]); }    
    inline T W(int q) const { return mW[q]; }
    inline int QRev(int q) const { return mQRev[q]; }
    inline T CSI() const { return mCSI; }
protected:
    static constexpr std::array<int, 6> mEX = {1, -1, 0, 0, 0, 0};
    static constexpr std::array<int, 6> mEY = {0, 0, 1, -1, 0, 0};
    static constexpr std::array<int, 6> mEZ = {0, 0, 0, 0, 1, -1};
    static constexpr std::array<T, 6> mW = {1./6., 1./6., 1./6., 1./6., 1./6., 1./6.};
    static constexpr std::array<int, 6> mQRev = {1, 0, 3, 2, 5, 4};
    const T mCSI = 2;
};

// Specialisation for D3Q7
template <typename T>
class AbstractLattice<T, 3, 7> : public AbstractLattice<T, 0, 0>
{
public:
    inline int EX(int q) const { return mEX[q]; }
    inline int EY(int q) const { return mEY[q]; }
    inline int EZ(int q) const { return mEZ[q]; }
    inline T CX(int q) const { return static_cast<T>(mEX[q]); }
    inline T CY(int q) const { return static_cast<T>(mEY[q]); }
    inline T CZ(int q) const { return static_cast<T>(mEZ[q]); }    
    inline T W(int q) const { return mW[q]; }
    inline int QRev(int q) const { return mQRev[q]; }
    inline T CSI() const { return mCSI; }
protected:
    static constexpr std::array<int, 7> mEX = {0, 1, -1, 0, 0, 0, 0};
    static constexpr std::array<int, 7> mEY = {0, 0, 0, 1, -1, 0, 0};
    static constexpr std::array<int, 7> mEZ = {0, 0, 0, 0, 0, 1, -1};
    static constexpr std::array<T, 7> mW = {1./4., 1./8., 1./8., 1./8., 1./8., 1./8., 1./8.};
    static constexpr std::array<int, 7> mQRev = {0, 2, 1, 4, 3, 6, 5};
    const T mCSI = 2;
};

// Specialisation for D3Q15
template <typename T>
class AbstractLattice<T, 3, 15> : public AbstractLattice<T, 0, 0>
{
public:
    inline int EX(int q) const { return mEX[q]; }
    inline int EY(int q) const { return mEY[q]; }
    inline int EZ(int q) const { return mEZ[q]; }
    inline T CX(int q) const { return static_cast<T>(mEX[q]); }
    inline T CY(int q) const { return static_cast<T>(mEY[q]); }
    inline T CZ(int q) const { return static_cast<T>(mEZ[q]); }    
    inline T W(int q) const { return mW[q]; }
    inline int QRev(int q) const { return mQRev[q]; }
    inline T CSI() const { return mCSI; }
protected:
    static constexpr std::array<int, 15> mEX = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1};
    static constexpr std::array<int, 15> mEY = {0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1};
    static constexpr std::array<int, 15> mEZ = {0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1};
    static constexpr std::array<T, 15> mW = {2./9., 1./9., 1./9., 1./9., 1./9., 1./9., 1./9., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72.};
    static constexpr std::array<int, 15> mQRev = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};
    const T mCSI = 3;
};

// Specialisation for D3Q15
template <typename T>
class AbstractLattice<T, 3, 19> : public AbstractLattice<T, 0, 0>
{
public:
    inline int EX(int q) const { return mEX[q]; }
    inline int EY(int q) const { return mEY[q]; }
    inline int EZ(int q) const { return mEZ[q]; }
    inline T CX(int q) const { return static_cast<T>(mEX[q]); }
    inline T CY(int q) const { return static_cast<T>(mEY[q]); }
    inline T CZ(int q) const { return static_cast<T>(mEZ[q]); }    
    inline T W(int q) const { return mW[q]; }
    inline int QRev(int q) const { return mQRev[q]; }
    inline T CSI() const { return mCSI; }
protected:
    static constexpr std::array<int, 19> mEX = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0};
    static constexpr std::array<int, 19> mEY = {0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1};
    static constexpr std::array<int, 19> mEZ = {0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1};
    static constexpr std::array<T, 19> mW = {1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,  1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.};
    static constexpr std::array<int, 19> mQRev = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};
    const T mCSI = 3;
};

// Specialisation for D3Q27
template <typename T>
class AbstractLattice<T, 3, 27> : public AbstractLattice<T, 0, 0>
{
public:
    inline int EX(int q) const { return mEX[q]; }
    inline int EY(int q) const { return mEY[q]; }
    inline int EZ(int q) const { return mEZ[q]; }
    inline T CX(int q) const { return static_cast<T>(mEX[q]); }
    inline T CY(int q) const { return static_cast<T>(mEY[q]); }
    inline T CZ(int q) const { return static_cast<T>(mEZ[q]); }    
    inline T W(int q) const { return mW[q]; }
    inline int QRev(int q) const { return mQRev[q]; }
    inline T CSI() const { return mCSI; }
protected:
    static constexpr std::array<int, 27> mEX = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1};
    static constexpr std::array<int, 27> mEY = {0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1};
    static constexpr std::array<int, 27> mEZ = {0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1};
    static constexpr std::array<T, 27> mW = {8./27., 2./27., 2./27., 2./27., 2./27., 2./27., 2./27., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216.};
    static constexpr std::array<int, 27> mQRev = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};
    const T mCSI = 3;
};
*/

#include "lattice/AbstractLattice.tpp"

#endif // ABSTRACTLATTICEDEF