/*****************************************************************************
 *                                                                           *
 *                            Struct definitions                             *
 *                                                                           *
 ****************************************************************************/

/**
 * @brief Characteristic scales for the Rayleigh Benard problem.
 * 
 * @param velocity The characteristic flow velocity scale.
 * @param time The characteristic flow development time scale.
 * @param length The characteristic flow length scale.
 */
struct CharacteristicScales
{
    const double velocity;
    const double time;
    const double length;

    // Constructor.
    CharacteristicScales(double u, double t, double l)
    : velocity(u), time(t), length(l) {}
};

/**
 * @brief Container for the dimensionless numbers governing the flow dynamics.
 * 
 * @param rayleigh The Rayleigh number.
 * @param prandtl The Prandtl number.
 */
struct DimensionlessNumbers 
{
    const double rayleigh;
    const double prandtl;

    // Constructor.
    DimensionlessNumbers(double ra, double pr) : rayleigh(ra), prandtl(pr) {}
};

/**
 * @brief Container for parameters related to the grid.
 * 
 * @param nx The number of grid points along the x-axis.
 * @param ny The number of grid points along the y-axis.
 * @param size The total number of grid points.
 */
struct GridParameters
{
    const int nx;
    const int ny;
    const int size;

    /**
     * @brief Constructor for GridParameters.
     * @param nx The number of grid points along the x-axis.
     * @param ny The number of grid points along the y-axis.
     */
    GridParameters(int nx, int ny) : nx(nx), ny(ny), size(nx * ny) {}
};

/**
 * @brief Parameters used in the simulation.
 * 
 * @param ag The product of thermal expansion coefficient (alpha)
 *        and gravitational acceleration (g).
 * @param omega_f The fluid relaxation rate in the LBGK collision operator.
 * @param omega_g The thermal relaxation rate in the LBGK collision operator.
 */
struct SimulationParameters 
{
    const double ag;
    const double omega_f;
    const double omega_g;

    // Constructor.
    SimulationParameters(double alpha_gravity, double relaxation_rate_fluid, double relxation_rate_thermal) 
        : ag(alpha_gravity), omega_f(relaxation_rate_fluid), omega_g(relxation_rate_thermal) {}
};