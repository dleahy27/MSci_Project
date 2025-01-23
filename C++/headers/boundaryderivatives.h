#pragma once
#include <cmath>
#include <vector>
#include <string>
#include "../headers/extrapolate.h"
#include "../headers/derivatives.h"

/// @brief Calculates the boundary derivatives needed in the bicubic spline algorithm.
///
/// This class calculates derivatives by implimenting the non-uniform finite difference alogorithm.
/// The higher order calculations use the same algorithm just using the previous orders grid values,
/// saving time on expanding everything out for an analytical term in terms of the grid z. This is done
/// as each specific derivative is needed for bicubic spline tests. As just stated, this is used in the bicubicspline
/// class to evaluate the derivatives at the edge of the grid, but currently also serves to compare analytical
/// derivatives to numerical for testing purposes
///
class BoundaryDerivatives{
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

    /// Calculate the second partial derivative with respect to y
    void partialDerivative2Y();
    /// Calculate the second partial derivative with respect to x
    void partialDerivative2X();

    /// Calculate the fourth partial derivative with respect to x and y
    void partialDerivative2X2Y();

    public:
    // Internal second derivative grid boundaries
    /// Second x partial derivatives on the left boundary
    std::vector<double> d2_left;
    /// Second x partial derivatives on the right boundary
    std::vector<double> d2_right;
    /// Second y partial derivatives on the bottom boundary
    std::vector<double> d2_bottom;
    /// Second y partial derivatives on the top boundary
    std::vector<double> d2_top;
    /// Corner derivatives
    std::vector<double> d4y2x2;

    /// Grid boundary derivatives struct
    Derivatives boundary_derivatives;

    /// Constructor
    BoundaryDerivatives(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<std::vector<double>>& Z);

    /// Compiles all of the boundary derivatives into the derivatives struct
    void SetBoundaryDerivatives();

    /// Outputs all the derivative grids as a csv
    void outputDerivs(std::string filename);
    
    void outputXDerivs(std::string filename);

    void outputYDerivs(std::string filename);

    void outputCornerDerivs(std::string filename);
};