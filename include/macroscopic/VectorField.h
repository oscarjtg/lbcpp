#ifndef VECTORFIELDDEF
#define VECTORFIELDDEF

#include <vector>
#include <string>

template <typename T>
class VectorField
{
public:
    VectorField(int size_x, int size_y, int size_z, std::string id_string = "u", std::string run_id = "test", std::string savepath = "output", int proc_num = 0);

    ~VectorField() = default;

    T GetValue(int a, int i, int j, int k) const;

    void SetValue(T value, int a, int i, int j, int k);

    T GetValueWrap(int a, int i, int j, int k) const;

    void SetValueWrap(T value, int a, int i, int j, int k);

    void AddToValue(T additive_value, int a, int i, int j, int k);

    void SetToConstantValue(T value, int a);

    void SetToConstantValue(T value);

    void SetLinearGradientZ(int a, T bottom, T top);

    T ComputeAverage(int a) const;

    T ComputeAverage() const;

    void DisplayInfo() const;

    void WriteToTextFile(const int timestep, int a) const;

    void WriteToTextFile(const int timestep = 0) const;

    void SetRunID(std::string run_id);

    inline int GetNX() const { return mSizeX; }
    inline int GetNY() const { return mSizeY; }
    inline int GetNZ() const { return mSizeZ; }

private:
    const int mSizeX, mSizeY, mSizeZ, mGridSize;

    const int mDims = 3;

    std::string mIDString, mRunID, mSaveDirectory;
    const int mProcessNumber;

    std::vector<T> mArray;

    inline int idx(int a, int i, int j, int k) const
    {
        return a + mDims * (i + mSizeX * (j + mSizeY * k));
        //return i + mSizeX * (j + mSizeY * (k + mSizeZ * a));
    }

    std::string construct_basepath(const int timestep, int a) const;
};

#include "macroscopic/VectorField.tpp"

#endif // VECTORFIELDDEF