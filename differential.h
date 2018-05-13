#ifndef DIFFERENTIAL_H
#define DIFFERENTIAL_H

/**
 * Store information on weights for interpolatory differential opertions (i.e. integral and derivative) on a fixed number of points
 *
 * Stores a mesh of [0,1] and weights wich can be used to compute \int_0^1 f(x) \de x and f'(y) for y a point of the mesh.
 */
class Differential
{
public:
    /**
     * Type of the mesh, affect mainly quadrature
     */
    enum class Type
    {
        Gauss,                                                                                                          //!< Gauss nodes
        ClenshawCurtis                                                                                                  //!< Chebyshev nodes
    };
    Differential(size_t npoints, Type type);                                                                            //!< initialize data for a given order
    ~Differential();                                                                                                    //!< initialize data for a given order
    double nodes(size_t index, double start=0, double end=1);                                                           //!< get index-th node in the mesh for the interval [start, end]
    double quadratureWeights(size_t index, double start=0, double end=1);                                               //!< get index-th quadrature weight, properly scaled for nodes in [start, end] 
    double differentiationWeights(size_t index, size_t point, double start=0, double end=1);                            //!< get index-th weight for approximating derivative in point-th point, properly scaled for nodes in [start, end]
    double evalPolynomial(size_t idx, double point, double start=0, double end=1);                                      //!< evaluato at point the idx-th lagrange polynomial for given nodes in [start,end]
    double* StealNodes();
    double* StealQuadratureWeights();
    double* StealDifferentiationWeights();
    const double* getNodes();
    const double* getQuadratureWeights();
    const double* getDifferentiationWeights();

private:
    size_t npoints;
    double *nodesx;
    double *qw;
    double *dw;
    double *w;
};

#endif // DIFFERENTIAL_H
