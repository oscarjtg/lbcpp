#ifndef MACROSCOPIC2DDEF
#define MACROSCOPIC2DDEF

template<class T> class Macroscopic2D
{
public:
    Macroscopic2D(int size_x, int size_y);
    ~Macroscopic2D();
    T GetDensity(int i, int j) const;
    T GetVelocityX(int i, int j) const;
    T GetVelocityY(int i, int j) const;
    T GetTemperature(int i, int j) const;
    T GetSalinity(int i, int j) const;
    void SetDensity(T r_, int i, int j);
    void SetVelocityX(T u_, int i, int j);
    void SetVelocityY(T v_, int i, int j);
    void SetTemperature(T t_, int i, int j);
    void SetSalinity(T s_, int i, int j);
private:
    const int mSizeX;
    const int mSizeY;
    const int mGridSize;
    T* mpR;
    T* mpU;
    T* mpV;
    T* mpT;
    T* mpS;
    inline int scalar_index(int i, int j) const
    {
        return i + mSizeX * j;
    }
};

#endif // MACROSCOPIC2DDEF