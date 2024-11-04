#ifndef EQUILIBRIADEF
#define EQUILIBRIADEF

/**
 * @brief Computes the first order equilibrium distribution function.
 *
 * @param density The fluid density at a given grid cell.
 * @param velocityProjection The projection of fluid velocity on a lattice vector.
 * @param latticeWeight The lattice weight corresponding to the fluid density.
 *
 * @return The correponding first order equilibrium distribution function value.
 */
template<class T>
inline T computeFirstOrderEquilibrium(
    T density, 
    T velocityProjection, 
    T latticeWeight)
{
    return density * latticeWeight * (static_cast<T>(1.0) + static_cast<T>(3.0) * velocityProjection);
}

/**
 * @brief Computes the second order equilibrium distribution function.
 *
 * @param density The fluid density at a given grid cell, calculated as 1 + sum df.
 * @param velocityProjection The projection of fluid velocity on a lattice vector.
 * @param velocityMagnitudeSquared The square of the fluid velocity.
 * @param latticeWeight The lattice weight corresponding to the fluid density.
 *
 * @return The correponding second order equilibrium distribution function value.
 */
template<class T>
inline T computeSecondOrderEquilibrium(
    T density, 
    T velocityProjection, 
    T velocityMagnitudeSquared,
    T latticeWeight)
{
    return density * latticeWeight * (static_cast<T>(1.0) + static_cast<T>(4.5) * velocityProjection * velocityProjection - static_cast<T>(1.5) * velocityMagnitudeSquared + static_cast<T>(3.0) * velocityProjection);
}

/**
 * @brief Computes the first order deviation distribution function equilibrium.
 *
 * @param density The fluid density at a given grid cell, calculated as 1 + sum df.
 * @param densityAnomaly The fluid density anomaly, sum df.
 * @param velocityProjection The projection of fluid velocity on a lattice vector.
 * @param latticeWeight The lattice weight corresponding to the fluid density.
 *
 * @return The correponding first order deviation distribution function equilibrium value.
 */
template<class T>
inline T computeFirstOrderDeviationEquilibrium(
    T density, 
    T densityAnomaly,
    T velocityProjection, 
    T latticeWeight)
{
    return density * latticeWeight * static_cast<T>(3.0) * velocityProjection + densityAnomaly * latticeWeight;
}

/**
 * @brief Computes the second order deviation distribution function equilibrium.
 *
 * @param density The fluid density at a given grid cell, calculated as 1 + sum df.
 * @param densityAnomaly The fluid density anomaly, sum df.
 * @param velocityProjection The projection of fluid velocity on a lattice vector.
 * @param velocityMagnitudeSquared The square of the fluid velocity.
 * @param latticeWeight The lattice weight corresponding to the fluid density.
 *
 * @return The correponding second order deviation distribution function equilibrium value.
 */
template<class T>
inline T computeSecondOrderDeviationEquilibrium(
    T density, 
    T densityAnomaly,
    T velocityProjection, 
    T velocityMagnitudeSquared,
    T latticeWeight)
{
    return density * latticeWeight * (static_cast<T>(4.5) * velocityProjection * velocityProjection - static_cast<T>(1.5) * velocityMagnitudeSquared + static_cast<T>(3.0) * velocityProjection) + densityAnomaly * latticeWeight;
}


#endif // EQUILIBRIADEF