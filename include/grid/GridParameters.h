#ifndef GRIDPARAMETERSDEF
#define GRIDPARAMETERSDEF

/**
 * @brief Container for parameters related to the grid.
 * 
 * @param nx The number of grid points along the x-axis.
 * @param ny The number of grid points along the y-axis.
 * @param size The total number of grid points.
 * @param isBottom Boolean that states whether subdomain includes bottom 
 *                 boundary of the full domain.
 * @param isTop Boolean that states whether subdomain includes top 
 *              boundary of the full domain.
 * @param isLeft Boolean that states whether subdomain includes left
 *               boundary of the full domain.
 * @param isRight Boolean that states whether subdomain includes right
 *                boundary of the full domain.
 */
struct GridParameters
{
    const int nx;
    const int ny;
    const int size;
    const int nSubDomsX;
    const int nSubDomsY;
    const int subDomCoordX;
    const int subDomCoordY;
    const bool isBottom;
    const bool isTop;
    const bool isLeft;
    const bool isRight;

    GridParameters(int nx, int ny, int nSubDomsX, int nSubDomsY, int subDomCoordX, int subDomCoordY, bool is_bottom, bool is_top, bool is_left, bool is_right) : 
    nx(nx), ny(ny), size(nx * ny), nSubDomsX(nSubDomsX), nSubDomsY(nSubDomsY), subDomCoordX(subDomCoordX), subDomCoordY(subDomCoordY), isBottom(is_bottom), isTop(is_top), isLeft(is_left), isRight(is_right) {}
};

#endif // GRIDPARAMETERSDEF