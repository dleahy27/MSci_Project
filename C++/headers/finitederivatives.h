#pragma once
#include <cmath>
#include <vector>
#include <string>
#include "../headers/extrapolate.h"
#include "../headers/derivatives.h"

/// @brief Calculates and stores a the partial derivatives of a grid up to the 3rd order.
///
/// This class calculates derivatives by implimenting the non-uniform finite difference alogorithm.
/// The higher order calculations use the same algorithm just using the previous orders grid values,
/// saving time on expanding everything out for an analytical term in terms of the grid z. This is done
/// as each specific derivative is needed for bicubic spline tests. As just stated, this is used in the bicubicspline
/// class to evaluate the derivatives at the edge of the grid, but currently also serves to compare analytical
/// derivatives to numerical for testing purposes
///
/// @todo Write method to only evaluate grid derivatives based off of grid extrapolation once testing is complete as otehr terms are uselsess in production
class FiniteDerivatives{
    private:
    /// Grid points in first dimension
    const std::vector<double> x;
    /// Grid points in second dimension
    const std::vector<double> y;
    /// Functional values at all xy grid points
    const std::vector<std::vector<double>> z;

    // grid dimensions
    /// Number of grid points in the first dimension
    const int m;
    /// Number of grid points in the second dimension
    const int n;

    /// Calculate all first dimensional derivatives
    void finiteD1();
    /// Calculate all second dimensional derivatives
    void finiteD2();
    /// Calculate all third dimensional derivatives
    void finiteD3();
    
    /// Calculate the first partial derivative with respect to x
    void partialDerivative1X();
    /// Calculate the first partial derivative with respect to y    
    void partialDerivative1Y();
    /// Calculate the second partial derivative with respect to y
    void partialDerivative2Y();
    /// Calculate the second partial derivative with respect to x
    void partialDerivative2X();
    /// Calculate the third partial derivative with respect to y
    void partialDerivative3Y();
    /// Calculate the third partial derivative with respect to x
    void partialDerivative3X();

    /// Calculate the fourth partial derivative with respect to x and y
    void partialDerivative2X2Y();

    public:
    // Internal derivative grids
    std::vector<std::vector<double>> d1x;
    std::vector<std::vector<double>> d2x;
    std::vector<std::vector<double>> d3x;

    std::vector<std::vector<double>> d1y;
    std::vector<std::vector<double>> d2y;
    std::vector<std::vector<double>> d3y;

    // Internal second derivative grid boundaries
    std::vector<double> d2_left;
    std::vector<double> d2_right;
    std::vector<double> d2_bottom;
    std::vector<double> d2_top;
    // Corner derivatives
    std::vector<double> d4y2x2;

    /// Grid boundary derivatives 
    Derivatives boundary_derivatives;

    /// Constructor
    FiniteDerivatives(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<std::vector<double>>& Z);

    /// Stores the left boundary derivatives from the grid
    void d2Left();

    /// Stores the right boundary derivatives from the grid
    void d2Right();

    /// Stores the top boundary derivatives from the grid
    void d2Top();

    /// Stores the bottom boundary derivatives from the grid
    void d2Bottom();

    /// Runs all the boundary storing functions
    void finiteBoundary();

    /// Compiles all of the boundary derivatives into the derivatives struct
    void SetBoundaryDerivatives();

    /// Outputs all the derivative grids as a csv
    void outputDerivs(std::string filename);
};