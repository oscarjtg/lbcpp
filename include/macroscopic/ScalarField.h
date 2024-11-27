#ifndef SCALARFIELDDEF
#define SCALARFIELDDEF

#include <vector>
#include <string>

template <typename T>
class ScalarField
{
public:
    ScalarField(int size_x, int size_y, int size_z, std::string id_string = "f", std::string run_id = "test", std::string savepath = "output", int proc_num = 0);

    ~ScalarField() = default;

    T GetValue(int i, int j, int k) const;

    void SetValue(T value, int i, int j, int k);

    T GetValueWrap(int i, int j, int k) const;

    void SetValueWrap(T value, int i, int j, int k);

    void AddToValue(T additive_value, int i, int j, int k);

    void SetToConstantValue(T value);

    void SetLinearGradientZ(T bottom, T top);

    T ComputeAverage() const;

    void DisplayInfo() const;

    void WriteToTextFile(const int timestep = 0) const;

    inline int GetNX() const { return mSizeX; }
    inline int GetNY() const { return mSizeY; }
    inline int GetNZ() const { return mSizeZ; }

private:
    const int mSizeX, mSizeY, mSizeZ, mGridSize;

    std::string mIDString, mRunID, mSaveDirectory;
    const int mProcessNumber;

    std::vector<T> mArray;

    inline int idx(int i, int j, int k) const
    {
        return i + mSizeX * j + mSizeX * mSizeY * k;
    }

    std::string construct_basepath(const int timestep) const;
};

#include "macroscopic/ScalarField.tpp"

#endif // SCALARFIELDDEF