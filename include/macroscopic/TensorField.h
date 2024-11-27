#ifndef TENSORFIELDDEF
#define TENSORFIELDDEF

#include <vector>
#include <string>

template <typename T>
class TensorField
{
public:
    TensorField(int size_x, int size_y, int size_z, std::string id_string = "u", std::string run_id = "test", std::string savepath = "output", int proc_num = 0);

    ~TensorField() = default;

    T GetValue(int a, int b, int i, int j, int k) const;

    void SetValue(T value, int a, int b, int i, int j, int k);

    T GetValueWrap(int a, int b, int i, int j, int k) const;

    void SetValueWrap(T value, int a, int b, int i, int j, int k);

    void AddToValue(T additive_value, int a, int b, int i, int j, int k);

    void SetToConstantValue(T value, int a, int b);

    void SetToConstantValue(T value);

    void SetLinearGradientZ(int a, int b, T bottom, T top);

    T ComputeAverage(int a, int b) const;

    T ComputeAverage() const;

    void DisplayInfo() const;

    void WriteToTextFile(const int timestep, int a, int b) const;

    void WriteToTextFile(const int timestep = 0) const;

    inline int GetNX() const { return mSizeX; }
    inline int GetNY() const { return mSizeY; }
    inline int GetNZ() const { return mSizeZ; }

private:
    const int mSizeX, mSizeY, mSizeZ, mGridSize;

    const int mRows = 3, mCols = 3;

    std::string mIDString, mRunID, mSaveDirectory;
    const int mProcessNumber;

    std::vector<T> mArray;

    inline int idx(int a, int b, int i, int j, int k) const
    {
        return a + mRows * (b + mCols * (i + mSizeX * (j + mSizeY * k)));
    }

    std::string construct_basepath(const int timestep, int a, int b) const;
};

#include "macroscopic/TensorField.tpp"

#endif // TensorFIELDDEF