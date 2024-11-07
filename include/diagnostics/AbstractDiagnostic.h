#ifndef ABSTRACTDIAGNOSTICDEF
#define ABSRTACTDIAGNOSITCDEF

template <class T>
class AbstractDiagnostic
{
public:
    AbstractDiagnostic() = default();

    ~AbstractDiagnostic() = default();

    T ComputeDiagnostic() = 0;
};

#endif // ABSTRACTDIAGNOSTICDEF