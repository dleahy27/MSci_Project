#pragma once
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include <stdexcept>
#include <typeinfo>
#include <iostream>
#include <fstream>


/// @brief Calculates and stores a 1D cubic spline, then evaluates at discrete value(s).
///
/// This class interpolates a 1D grid of data by employing the standard cubic spline algorithm,
/// splitting up the domain into smaller sub domains that each have an interpolated cubic polynomial associated with them.
/// The coefficients of each cubic spline for each sub-domain are stored internally in a 2D vector and subsequent
/// points within the grid's domain can then be used to approximate the functional value in between grid points. 
class CubicSpline{
    private:
    /// Array of grid points
    std::vector<double> ts;
    /// Array of function values at grid points
    std::vector<double> ys;

    /// Number of grid points
    int N;

    /// Second derivative of the first/start grid point
    double d2t2_start;
    /// Second derivative of the last/end grid point 
    double d2t2_end;

    /// Internal 2D array containing spline coefficients at each interval
    std::vector<std::vector<double>> _S;

    /// Evaluates the functional value for a point within the grid
    double piecewise(const double t, const int it, std::function<double(double,  int)>& lambda);

    /// Evaluates an array of functional values for an array of points within the grid
    std::vector<double> piecewise(const std::vector<double>& t, const std::vector<int>& its, std::function<double(double,  int)>& lambda);

    public:
    /// Constructor
    CubicSpline(const std::vector<double>& Ts, const std::vector<double>& Ys, double D2t2_start, double D2t2_end);

    /// Evaluates a functional value at a specific point within the grid
    double evaluateSpline(double t, const std::vector<std::vector<double>>& S, const std::vector<double>& ts);

    /// Evaluates a vector of functional values for a vector of points within the grid
    std::vector<double> evaluateSpline(const std::vector<double>& t, const std::vector<std::vector<double>>& S, const std::vector<double>& ts);

    /// Get 2D spline vector
    std::vector<std::vector<double>> GetS();

    /// Set 2D spline vector
    void SetS(const std::vector<std::vector<double>>& S);
    
};

