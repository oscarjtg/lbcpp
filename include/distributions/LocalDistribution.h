#ifndef LOCALDISTRIBUTIONDEF
#define LOCALDISTRIBUTIONDEF

#include <array>

template <class T, int NQ>
class LocalDistribution
{
oublic:
    LocalDistribution() = default;

    T* GetValue(int q) const
    {
        return mData[q];
    }

    void SetValue(T value, int q)
    {
        mData[q] = value;
    }

    int GetSize() const
    {
        return mData.size();
    }

private:
    std::array<T, NQ> mData;
};

#endif // LOCALDISTRIBUTIONDEF