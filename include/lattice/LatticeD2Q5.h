#ifndef LATTICED2Q5DEF
#define LATTICED2Q5DEF

template<class T> class LatticeD2Q5
{
public:
    LatticeD2Q5(int size_x, int size_y);
    ~LatticeD2Q5();
    T GetF0(int i, int j) const;
    T GetF1(int i, int j) const;
    T GetF2(int i, int j) const;
    T GetF3(int i, int j) const;
    T GetF4(int i, int j) const;
    T GetPreStreamF0(int i, int j) const;
    T GetPreStreamF1(int i, int j) const;
    T GetPreStreamF2(int i, int j) const;
    T GetPreStreamF3(int i, int j) const;
    T GetPreStreamF4(int i, int j) const;
    T GetBouncedF0(int i, int j) const;
    T GetBouncedF1(int i, int j) const;
    T GetBouncedF2(int i, int j) const;
    T GetBouncedF3(int i, int j) const;
    T GetBouncedF4(int i, int j) const;
    void SetF0(T f_, int i, int j);
    void SetF1(T f_, int i, int j);
    void SetF2(T f_, int i, int j);
    void SetF3(T f_, int i, int j);
    void SetF4(T f_, int i, int j);
    void SetInitF0(T f_, int i, int j);
    void SetInitF1(T f_, int i, int j);
    void SetInitF2(T f_, int i, int j);
    void SetInitF3(T f_, int i, int j);
    void SetInitF4(T f_, int i, int j);
    void SwapPointers();
    void DisplayPointerInfo() const;
    void DisplayInfo() const;
    void WriteToCSV(const std::string& path, const char letter, const std::string& runId = "test", const int timestep = 0) const;
private:
    const int mSizeX;
    const int mSizeY;
    const int mGridSize;
    T* mpF0;
    T* mpF1;
    T* mpF2;
    T* mpF3;
    T* mpF4;
    inline int scalar_index(int i, int j) const
    {
        return i + mSizeX * j;
    }
    std::string construct_basepath(const std::string& path, const std::string& runId, const int timestep, const char letter) const;
    void write_csv(const std::string& path, T (LatticeD2Q5<T>::*get_func)(int, int) const) const;
};

#endif // LATTICED2Q5DEF