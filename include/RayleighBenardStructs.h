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
    const float velocity;
    const float time;
    const float length;

    // Constructor.
    CharacteristicScales(float u, float t, float l)
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
    const float rayleigh;
    const float prandtl;

    // Constructor.
    DimensionlessNumbers(float ra, float pr) : rayleigh(ra), prandtl(pr) {}
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
    const float ag;
    const float omega_f;
    const float omega_g;

    // Constructor.
    SimulationParameters(float alpha_gravity, float relaxation_rate_fluid, float relxation_rate_thermal) 
        : ag(alpha_gravity), omega_f(relaxation_rate_fluid), omega_g(relxation_rate_thermal) {}
};