#pragma once
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include "../headers/derivatives.h"

namespace funcs{
    std::vector<double> random_uniform(double start, double end, unsigned int N);

    double gauss(const double& x, const double& y);

    std::vector<std::vector<double>> gauss(const std::vector<double>& x, const std::vector<double>& y);

    std::vector<std::vector<double>> trig(const std::vector<double>& x, const std::vector<double>& y);

    std::vector<std::vector<double>> pdfLike(const std::vector<double>& x, const std::vector<double>& y);

    std::vector<double> d2dx2_gauss(const double x, const std::vector<double>& y);

    std::vector<double> d2dy2_gauss(const std::vector<double>& x, const double y);

    std::vector<double> d4dx2dy2_gauss(const std::vector<double>& X, const std::vector<double>& Y);

    Derivatives deriv_gauss(const std::vector<double>& x, const std::vector<double>& y);

    std::vector<double> d1dx_trig(const double x, const std::vector<double>& y);

    std::vector<double> d1dy_trig(const std::vector<double>& x, const double y);

    std::vector<double> d2dx2_trig(const double x, const std::vector<double>& y);

    std::vector<double> d2dy2_trig(const std::vector<double>& x, const double y);

    std::vector<double> d4dx2dy2_trig(const std::vector<double>& X, const std::vector<double>& Y);

    Derivatives deriv_trig(const std::vector<double>& x, const std::vector<double>& y);

    std::vector<std::vector<double>> d1dx_pdfLike(const std::vector<double>& x, const std::vector<double>& y);

    std::vector<std::vector<double>> d1dy_pdfLike(const std::vector<double>& x, const std::vector<double>& y);

    std::vector<std::vector<double>> d2dx2_pdfLike(const std::vector<double>& x, const std::vector<double>& y);

    std::vector<std::vector<double>> d2dy2_pdfLike(const std::vector<double>& x, const std::vector<double>& y);

    std::vector<double> d4dx2dy2_pdfLike(const std::vector<double>& X, const std::vector<double>& Y);

    Derivatives deriv_pdfLike(const std::vector<double>& x, const std::vector<double>& y);

    void outputTrigDerivatives( const std::vector<double>& x, const std::vector<double>& y, const std::string& filename );

    void outputPdfDerivatives( const std::vector<double>& x, const std::vector<double>& y, const std::string& filename);
}