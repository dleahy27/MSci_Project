#pragma once
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include <stdexcept>
#include <typeinfo>
#include <iostream>
#include <fstream>


class CubicSpline{
    private:
    std::vector<double> ts;
    std::vector<double> ys;

    int N;

    double d2t2_start;
    double d2t2_end;

    std::vector<std::vector<double>> _S;

    double piecewise(double x, std::vector<bool> conds, std::vector<std::function<double(double,  int)>> lambda);

    std::vector<double> piecewise(std::vector<double> x, std::vector<std::vector<bool>> conds, std::vector<std::function<double(double,  int)>> lambda);

    public:
    // Constructor
    CubicSpline(const std::vector<double>& Ts, const std::vector<double>& Ys, const double& D2t2_start, const double& D2t2_end);

    // Spline evaluation at value(s) x
    double evaluateSpline(double x,std::vector<std::vector<double>>S, std::vector<double>ts, bool verbose);

    std::vector<double> evaluateSpline(std::vector<double> x,std::vector<std::vector<double>>S, std::vector<double>ts, bool verbose);

    // Getters
    std::vector<std::vector<double>> GetS();
    int GetN();

    // Setters
    void SetS(std::vector<std::vector<double>> S);
    
};

