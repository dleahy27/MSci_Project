#pragma once
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <Eigen/Dense>
#include "../headers/cubicspline.h"
#include "../headers/derivatives.h"

// define matrix classes from eigen
/// Eigen matrix class for a 16x16 array
typedef Eigen::Matrix<double, 16, 16> Matrix16;
/// Eigen vector class for a size 16 array
typedef Eigen::Vector<double, 16> Vector16;

/// @brief Calculates and stores a 2D bicubic spline, then evaluates at discrete value(s).
///
/// This class interpolates a 2D grid of data by employing the standard bicubic spline algorithm,
/// splitting up the domain into smaller 2D sub-domains that each have an interpolated bicubic polynomial associated with them.
/// The coefficients of each bicubic cubic spline for each sub-domain are stored internally in a 3D vector and subsequent
/// points within the grid's domain can then be used to approximate the functional value in between grid points. 
class BicubicSpline{
    private:
    // variables

    // parameter vectors
    // Grid interpolated in logspace i.e. u = log x and v = log Q^2
    /// Grid points in first dimension
    std::vector<double>us;
    /// Grid points in second dimension
    std::vector<double>vs;
    /// Functional values at the grid points
    std::vector<std::vector<double>>xfs;

    // x-axis points to evaluate derivatives
    std::vector<double>dus;
    /// v-axis points to evaluate derivatives
    std::vector<double>dvs;
    // Derivative vectors -- differentiability and continuity checks
    std::vector<std::vector<double>>dxfdx_vec;
    std::vector<std::vector<double>>d2xfd2x_vec;
    std::vector<std::vector<double>>d3xfd3x_vec;

    std::vector<std::vector<double>>dxfdQ2_vec;
    std::vector<std::vector<double>>d2xfd2Q2_vec;
    std::vector<std::vector<double>>d3xfd3Q2_vec;

    // grid dimensions
    /// Number of first dimension grid points
    unsigned int m;
    /// Number of second dimension grid points
    unsigned int n;

    // Derivatives grid dimensions
    unsigned int dm;
    unsigned int dn;

    // second derivative arrays
    /// Second partial x derivative at the left boundary
    std::vector<double> d2xfd2xs_left;
    /// Second partial x derivative at the right boundary
    std::vector<double> d2xfd2xs_right;
    /// Second partial Q^2 derivative at the bottoom boundary
    std::vector<double> d2xfd2Q2s_bottom;
    /// Second partial Q^2 derivative at the top boundary
    std::vector<double> d2xfd2Q2s_top;
    /// Fourth partial mix derivative at the corners
    std::vector<double> d4xfd2xd2Q2s_corners;

    /// Internal spline coefficients array
    std::vector<double> S;

    // functions
    /// Calculates the bicubic spline coefficients for the grid
    void CalculateSpline();

    // Below calculates derivatives, first are the depreciated methods used for the project
    // and after that are the actual methods that are to be used in LHAPDF
    /// Initialise output derivative arrays -- used for project
    void initDerivs(const std::vector<double>& u, const std::vector<double>& v);

    /// Derivative calculations -- used for project

    void dxfdx(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<int>& ius, const std::vector<int>& ivs);

    void d2xfd2x(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<int>& ius, const std::vector<int>& ivs);

    void d3xfd3x(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<int>& ius, const std::vector<int>& ivs);

    void dxfdQ2(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<int>& ius, const std::vector<int>& ivs);

    void d2xfd2Q2(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<int>& ius, const std::vector<int>& ivs);

    void d3xfd3Q2(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<int>& ius, const std::vector<int>& ivs);

    /// Derivative evaluations -- for LHAPDF

    double dxfdx(double u, double v, size_t ix, size_t iy) const;

    double d2xfd2x(double u, double v, size_t ix, size_t iy) const;

    double dxfdQ2(double u, double v, size_t ix, size_t iy) const;

    double d2xfd2Q2(double u, double v, size_t ix, size_t iy) const;

    double dxfdu(double u, double v, size_t ix, size_t iy) const;

    double d2xfd2u(double u, double v, size_t ix, size_t iy) const;

    double dxfdv(double u, double v, size_t ix, size_t iy) const;

    double d2xfd2v(double u, double v, size_t ix, size_t iy) const;

    public:

    /// Default constructor, for storage in STL containers
    // BicubicSpline() {}

    /// Standard constructor
    BicubicSpline(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<std::vector<double>>& xfs);

    /// Constructor with the additional argument providing boundary derivatives
    BicubicSpline(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<std::vector<double>>& xfs, const Derivatives& ds);

    /// Evaluates the bicubic spline at a singular xy point within the grid
    double evaluateSpline( double u, double v ) const;

    /// Evaluates the bicubic spline at a singular xy point within the grid when the grid indices are also known
    double evaluateSpline( double u, double v, size_t iu, size_t iv ) const;
    
    /// Evaluates the bicubic spline at for an array of xy point within the grid
    std::vector<std::vector<double>> evaluateSpline( const std::vector<double>& u, const std::vector<double>& v ) const;

    /// Depreciated method to calculate array of derivative values -- used in project
    void calculateDerivs(const std::vector<double>& u, const std::vector<double>& v);

    /// LHAPDF method to evaluate the x-Q2 derivatives at a log grid point 
    std::tuple<double,double,double,double> evaluateDerivs(double du, double dy, size_t idu, size_t idv) const;
    
    /// LHAPDF method to evaluate the u-v derivatives at a log grid point 
    std::tuple<double,double,double,double> evaluateLogDerivs(double du, double dy, size_t idu, size_t idv) const;

    /// Depreciated method to write the derivative vectors out to file -- used in project
    void outputDerivs(const std::string& filename) const;

};
