#ifndef ABSTRACTDIAGNOSTICDEF
#define ABSRTACTDIAGNOSITCDEF

template <typename T>
class AbstractDiagnostic
{
public:
    AbstractDiagnostic() = default;

    ~AbstractDiagnostic() = default;

    double ComputeDiagnostic() = 0;
};

#endif // ABSTRACTDIAGNOSTICDEF