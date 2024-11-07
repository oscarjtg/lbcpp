#ifndef FIXEDSIZECONTAINERDEF
#define FIXEDSIZECONTAINERDEF

#include <iostream>
#include <deque>
#include <numeric> // for accumulate

template <class T>
class FixedSizeContainer {
public:
    FixedSizeContainer(int size) : mSize(size) {}

    ~FixedSizeContainer() = default;

    void AddElement(T element) 
    {
        if (static_cast<int>(mContainer.size()) == mSize) 
        {
            mContainer.pop_back(); // Remove the last element if the container is full
        }
        mContainer.push_front(element); // Add the new element to the front
    }

    T GetElement(int i)
    {
        if (i < static_cast<int>(mContainer.size()))
        {
            return mContainer[i];
        }
        else
        {
            return 0;
        }
    }

    void PrintContainer() const 
    {
        for (const auto& elem : mContainer) 
        {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    T CalculateAverage() const {
        if (mContainer.empty())
        {
            return 0.0;
        }
        double sum = std::accumulate(mContainer.begin(), mContainer.end(), 0.0);
        return sum / mContainer.size();
    }

private:
    std::deque<T> mContainer;
    int mSize;
};

#endif // FIXEDSIZECONTAINERDEF