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

    double piecewise(const double t, const int it, std::function<double(double,  int)>& lambda);

    std::vector<double> piecewise(const std::vector<double>& t, const std::vector<int>& its, std::function<double(double,  int)>& lambda);

    public:
    // Constructor
    CubicSpline(const std::vector<double>& Ts, const std::vector<double>& Ys, double D2t2_start, double D2t2_end);

    // Spline evaluation at value(s) x
    double evaluateSpline(double t, const std::vector<std::vector<double>>& S, const std::vector<double>& ts);

    std::vector<double> evaluateSpline(const std::vector<double>& t, const std::vector<std::vector<double>>& S, const std::vector<double>& ts);

    // Getters
    std::vector<std::vector<double>> GetS();
    
    int GetN();

    // Setters
    void SetS(const std::vector<std::vector<double>>& S);
    
};

