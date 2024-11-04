#ifndef MACROSCOPICVARIABLEDEF
#define MACROSCOPICVARIABLEDEF

#include <vector>
#include <string>

template <class T>
class MacroscopicVariable
{
public:
    MacroscopicVariable(int size_x, int size_y, int size_z, std::string id_string = "f", std::string run_id = "test", std::string savepath = "output", int proc_num = 0);

    ~MacroscopicVariable() = default;

    T GetValue(int i, int j, int k) const;

    void SetValue(T value, int i, int j, int k);

    void SetLinearGradientZ(T bottom, T top);

    T ComputeAverage() const;

    void DisplayInfo() const;

    void WriteToTextFile(const int timestep = 0) const;

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

#include "macroscopic/MacroscopicVariable.tpp"

#endif // MACROSCOPICVARIABLEDEF