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
/// saving time on expanding everything out for an analytical term in terms of the grid xf. This is done
/// as each specific derivative is needed for bicubic spline tests. As just stated, this is used in the bicubicspline
/// class to evaluate the derivatives at the edge of the grid, but currently also serves to compare analytical
/// derivatives to numerical for testing purposes
///
class BoundaryDerivatives{
    private:
    /// Grid points in first dimension
    const std::vector<double> us;
    /// Grid points in second dimension
    const std::vector<double> vs;
    /// Functional values at all xy grid points
    const std::vector<std::vector<double>> xfs;

    // grid dimensions
    /// Number of grid points in the first dimension
    const int m;

    /// Number of grid points in the second dimension
    const int n;

    /// Calculate the second partial derivative with respect to Q2
    void d2xfd2Q2();

    /// Calculate the second partial derivative with respect to x
    void d2xfd2x();

    /// Calculate the fourth partial derivative with respect to x and Q2
    void d4xfd2xd2Q2();

    public:
    // Internal second derivative grid boundaries
    /// Second x partial derivatives on the left boundary
    std::vector<double> d2xfd2x_left;
    /// Second x partial derivatives on the right boundary
    std::vector<double> d2xfd2x_right;
    /// Second Q2 partial derivatives on the bottom boundary
    std::vector<double> d2xfd2Q2_bottom;
    /// Second Q2 partial derivatives on the top boundary
    std::vector<double> d2xfd2Q2_top;
    /// Corner derivatives
    std::vector<double> d4xfd2xd2Q2_vec;

    // Internal grid values
    /// Left boundary xf
    std::vector<double> xfs_left;
    /// Right boundary xf
    std::vector<double> xfs_right;
    /// Top boundary xf
    std::vector<double> xfs_top;
    /// Bottom boundary xf
    std::vector<double> xfs_bottom;

    /// Grid boundary derivatives struct
    Derivatives boundary_derivatives;

    /// Constructor
    BoundaryDerivatives(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<std::vector<double>>& xfs);

    /// Compiles all of the boundary derivatives into the derivatives struct
    void SetBoundaryDerivatives();

    /// Outputs all the derivative grids as a csv
    void outputDerivs(std::string filename);
    
    /// Outputs x-derivative grids as a csv
    void outputXDerivs(std::string filename);

    /// Outputs Q2-derivative grids as a csv
    void outputQ2Derivs(std::string filename);

    /// Outputs corner-derivative grids as a csv
    void outputCornerDerivs(std::string filename);

    /// Outputs all the extrapolated xf boundaries as a csv
    void outputxfs(std::string filename);
    
    /// Outputs the left and right extrapolated xf boundaries as a csv
    void outputLeftRightxfs(std::string filename);

    /// Outputs the top and bottom extrapolated xf boundaries as a csv
    void outputTopBottomxfs(std::string filename);
};