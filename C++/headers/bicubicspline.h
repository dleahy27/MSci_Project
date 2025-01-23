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
    /// Grid points in first dimension
    std::vector<double>xs;
    /// Grid points in second dimension
    std::vector<double>ys;
    /// Functional values at the grid points
    std::vector<std::vector<double>>zs;

    // x-axis points to evaluate derivatives
    std::vector<double>dxs;
    /// y-axis points to evaluate derivatives
    std::vector<double>dys;
    // Derivative vectors -- differentiability and continuity checks
    std::vector<std::vector<double>>d1x;
    std::vector<std::vector<double>>d2x;
    std::vector<std::vector<double>>d3x;

    std::vector<std::vector<double>>d1y;
    std::vector<std::vector<double>>d2y;
    std::vector<std::vector<double>>d3y;

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
    std::vector<double> d2x2s_left;
    /// Second partial x derivative at the right boundary
    std::vector<double> d2x2s_right;
    /// Second partial y derivative at the bottoom boundary
    std::vector<double> d2y2s_bottom;
    /// Second partial y derivative at the top boundary
    std::vector<double> d2y2s_top;
    /// Fourth partial mix derivative at the corners
    std::vector<double> d4x2y2s_corners;

    /// Internal spline coefficients array
    std::vector<std::vector<Vector16>> S;

    // functions
    /// Calculates the bicubic spline coefficients for the grid
    void CalculateSpline();

    void initDerivs(const std::vector<double>& X, const std::vector<double>& Y);

    // Derivative calculations
    void d1X(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<int>& iXs, const std::vector<int>& iYs);

    void d2X(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<int>& iXs, const std::vector<int>& iYs);

    void d3X(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<int>& iXs, const std::vector<int>& iYs);

    void d1Y(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<int>& iXs, const std::vector<int>& iYs);

    void d2Y(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<int>& iXs, const std::vector<int>& iYs);

    void d3Y(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<int>& iXs, const std::vector<int>& iYs);

    public:
    /// Standard constructor
    BicubicSpline(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<std::vector<double>>& Zs);

    /// Constructor with the additional argument providing boundary derivatives
    BicubicSpline(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<std::vector<double>>& Zs, const Derivatives& ds);

    /// Evaluates the bicubic spline at a singular xy point within the grid
    double evaluateSpline( double X, double Y );
    /// Evaluates the bicubic spline at for an array of xy point within the grid
    std::vector<std::vector<double>> evaluateSpline( const std::vector<double>& X, const std::vector<double>& Y );

    void calculateDerivs(const std::vector<double>& X, const std::vector<double>& Y);

    void outputDerivs(const std::string& filename);
};