#pragma once
#include "linalg.h"
#include "cubicspline.h"
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <algorithm>

class BicubicSpline{
    private:
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
    std::vector<std::vector<std::vector<double>>> _S;

    public:
    // eventually split this constructor into seperate methods
    BicubicSpline(std::vector<double>Xs, std::vector<double>Ys, std::vector<std::vector<double>>Zs);

    double evaluateSpline( double X, double Y );

    std::vector<std::vector<double>> evaluateSpline( std::vector<double> X, std::vector<double> Y );

    void SetS(std::vector<std::vector<std::vector<double>>> S);

    std::vector<std::vector<std::vector<double>>> GetS();
};