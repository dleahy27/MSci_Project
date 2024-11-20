#pragma once
#include <cmath>
#include <vector>
#include <string>
#include "../headers/extrapolate.h"
#include "../headers/derivatives.h"

class FiniteDerivatives{
    private:
    const std::vector<double> x;
    const std::vector<double> y;

    const std::vector<std::vector<double>> z;

    const int m;
    const int n;

    void finiteD1();
    
    void finiteD2();
    
    void partialDerivative1X();

    void partialDerivative1Y();

    void partialDerivative2Y();

    void partialDerivative2X();

    void partialDerivative2X2Y();

    public:
    std::vector<std::vector<double>> d1x;
    std::vector<std::vector<double>> d2x;

    std::vector<std::vector<double>> d1y;
    std::vector<std::vector<double>> d2y;

    std::vector<double> d2_left;
    std::vector<double> d2_right;
    std::vector<double> d2_bottom;
    std::vector<double> d2_top;

    std::vector<double> d4y2x2;

    Derivatives boundary_derivatives;

    FiniteDerivatives(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<std::vector<double>>& Z);

    void d2Left();

    void d2Right();

    void d2Top();

    void d2Bottom();

    void finiteBoundary();

    void SetBoundaryDerivatives();

    void outputDerivs(std::string filename);
};