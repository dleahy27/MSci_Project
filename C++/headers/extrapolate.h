#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>

// struct
struct Coeffs{
    double x0;
    double y0;
    double a;
    double b;
};

/// @brief Extrapolates the functional values of a 1D grid of data beyond the grid.
///
/// This class fits the 1D grid data for the 2 or 3 point slice provided and then extrapolates to a point beyond the grid.
/// Has two current extrapolation schemes: a linear model (2 points) and a power law model (3 points).
/// Used within the finitederivatives class to extend the second order partial derivatives to the boundary values,
/// the current finite difference scheme employed cannot evaluate these.
///
/// @todo Add a log_10 fit and evaluation for better fitting along specific boundaries

class Extrapolate{
    private:
    /// Array of grid points
    const std::vector<double> xs;
    /// Array of function values at grid points
    const std::vector<double> ys;

    /// Evaluates the second derivative of the fitted power law at a new point
    double secondDerivPower(double x, const Coeffs& coeffs);

    /// Evaluates the second derivative of the fitted ten power law at a new point
    double secondDerivTenPower(double x, const Coeffs& coeffs);

    /// Evaluates the value of the power law model at a point
    double powerLawModel(double x, const Coeffs& coeffs);

    /// Evaluates the value of the ten power law model at a point
    double tenPowerLawModel(double x, const Coeffs& coeffs);

    /// Evaluates the value of the linear model at a point
    double linearModel(double x, const Coeffs& coeffs);

    /// Fits the power law model based off of grid points and function values
    Coeffs interpPowerLaw();
    /// Fits the 10 power law model based off of grid points and function values
    Coeffs interpTenPowerLaw();
    /// Fits the linear model based off of grid points and function values
    Coeffs interpLinearLaw();

    

    public:
    /// Extrapolated value
    double solution;

    /// Constructor
    Extrapolate(const std::vector<double>& X, const std::vector<double>& Y);

    /// Power law extrapolation beyond the grid for point x 
    void extrapolatePowerLaw(double x);
    /// 10 Power law extrapolation beyond the grid for point x 
    void extrapolateTenPowerLaw(double x);
    /// Linear extrapolation beyond the grid for point x 
    void extrapolateLinearLaw(double x);

    /// Output the extrapolated z values (before taking the second deriv)
    double tenPowerLawZs(double x);

    double powerLawZs(double x);

};