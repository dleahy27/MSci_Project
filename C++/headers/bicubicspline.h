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

typedef Eigen::Matrix<double, 16, 16> Matrix16;
typedef Eigen::Vector<double, 16> Vector16;

class BicubicSpline{
    private:
    // variables

    // parameter vectors
    std::vector<double>xs;
    std::vector<double>ys;
    std::vector<std::vector<double>>zs;

    // grid dimensions
    int m;
    int n;

    // second derivative arrays
    std::vector<double> d2x2s_left;
    std::vector<double> d2x2s_right;
    std::vector<double> d2y2s_bottom;
    std::vector<double> d2y2s_top;
    std::vector<double> d4x2y2s_corners;

    // internal spline grids vector
    std::vector<std::vector<Vector16>> S;

    // functions

    void CalculateSpline();

    public:
    // Constructors
    BicubicSpline(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<std::vector<double>>& Zs);

    BicubicSpline(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<std::vector<double>>& Zs, const Derivatives& ds);

    // Functions
    double evaluateSpline( double X, double Y );

    std::vector<std::vector<double>> evaluateSpline( const std::vector<double>& X, const std::vector<double>& Y );
};