#ifndef LATTICED2Q9DEF
#define LATTICED2Q9DEF

template<class T> class LatticeD2Q9
{
public:
    LatticeD2Q9(int size_x, int size_y);
    ~LatticeD2Q9();
    T GetF0(int i, int j) const;
    T GetF1(int i, int j) const;
    T GetF2(int i, int j) const;
    T GetF3(int i, int j) const;
    T GetF4(int i, int j) const;
    T GetF5(int i, int j) const;
    T GetF6(int i, int j) const;
    T GetF7(int i, int j) const;
    T GetF8(int i, int j) const;
    void SetF0(T f_, int i, int j);
    void SetF1(T f_, int i, int j);
    void SetF2(T f_, int i, int j);
    void SetF3(T f_, int i, int j);
    void SetF4(T f_, int i, int j);
    void SetF5(T f_, int i, int j);
    void SetF6(T f_, int i, int j);
    void SetF7(T f_, int i, int j);
    void SetF8(T f_, int i, int j);
    void SetF9(T f_, int i, int j);
    void SetInitF0(T f_, int i, int j);
    void SetInitF1(T f_, int i, int j);
    void SetInitF2(T f_, int i, int j);
    void SetInitF3(T f_, int i, int j);
    void SetInitF4(T f_, int i, int j);
    void SetInitF5(T f_, int i, int j);
    void SetInitF6(T f_, int i, int j);
    void SetInitF7(T f_, int i, int j);
    void SetInitF8(T f_, int i, int j);
    void SwapPointers();
    void DisplayPointerInfo() const;
    void DisplayInfo() const;
private:
    const int mSizeX;
    const int mSizeY;
    const int mGridSize;
    T* mpF0;
    T* mpF1;
    T* mpF2;
    T* mpF3;
    T* mpF4;
    T* mpF5;
    T* mpF6;
    T* mpF7;
    T* mpF8;
    inline int scalar_index(int i, int j) const
    {
        return i + mSizeX * j;
    }
};

#endif // LATTICED2Q9DEF