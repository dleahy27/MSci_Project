#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

// struct
struct Coeffs{
    double x0;
    double y0;
    double a;
    double b;
};

class Extrapolate{
    private:
    const std::vector<double> xs;
    const std::vector<double> ys;

    std::vector<double> gradients;

    double powerLawModel(double x, const Coeffs& coeffs);

    Coeffs interpPowerLaw();
    Coeffs interpLinearLaw();

    double extrapolatePowerLaw(double x);
    double extrapolateLinearLaw(double x);

    public:

    //variables
    double solution;

    // constructor
    Extrapolate(double x, const std::vector<double>& X, const std::vector<double>& Y);

    double linearModel(double x, const Coeffs& coeffs);

};