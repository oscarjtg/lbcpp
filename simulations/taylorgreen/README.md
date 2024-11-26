# Taylor Green Vortex Decay Benchmark

This study tests the decay of the Taylor Green Vortex over one decay timescale.

Compare different collision operators (currently SRT, TRT)
and initialisation schemes (feq, feq+fneq, Wei's consistent initialisation scheme)
and floating point precision (F64, F32).

Output files contain columns for timestep, error in r, error in u, and error in v.

The error is computed as the L2 error norm between the computed fields and the analytic solution.

The files are named as such:

{runid}_l2error.csv

{runid} is the ID of the run. It takes the following form:

{runid} = {collision_operator}_{initialisation_scheme}_{dfType}

{collision_operator}:
    SRTxx : single relaxation time
    TRT04 : two relaxation time, with magic parameter 1/4
    TRT12 : two relaxation time, with magic parameter 1/12

{initialisation_scheme}:
    FEQ : equilibrium initialisation scheme
    NEQ : non-equilibrium initialisation scheme
    WEI : Wei's consistent initialisation scheme

{dfType}:
    F64 : double precision
    F32 : single precision