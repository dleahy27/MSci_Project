#include "../headers/extrapolate.h"

// fast square function
inline double square(double x){return x*x;}

// Model function y = a(x - x0)^b + y0
double Extrapolate::powerLawModel(double x, const Coeffs& coeffs) {
    // Covers some backwards cases by reflecting around x-x0
    if ( (x - coeffs.x0)<0 ){
        return (coeffs.y0 + coeffs.a*std::pow((-x + coeffs.x0),coeffs.b));
    }else{
        return (coeffs.y0 + coeffs.a*std::pow(((x - coeffs.x0)),coeffs.b));
    }
    
}

// Analytical second derivative of power function
double Extrapolate::secondDerivPower(double x, const Coeffs& coeffs) {
    // Covers some backwards cases by reflecting around x-x0
    if (coeffs.b == 1){
        return 0;
    } else if ( x - coeffs.x0<0 ){
        return -coeffs.a*coeffs.b * ( -(coeffs.b - 1)*std::pow(-x + coeffs.x0, coeffs.b - 2) - std::log(10)*std::pow((coeffs.x0 - x),coeffs.b - 1) ) / (std::pow(10,2*x)*std::log(10)*std::log(10));
    }else if ( x - coeffs.x0>0 ){
        return  coeffs.a*coeffs.b * ( (coeffs.b - 1)*std::pow(x - coeffs.x0, coeffs.b - 2) - std::log(10)*std::pow((x - coeffs.x0),coeffs.b - 1) ) / (std::pow(10,2*x)*std::log(10)*std::log(10));
    }else {
        return 0;
    }
}

// Linear model: y = a(x-x0) + y0
double Extrapolate::linearModel(double x, const Coeffs& coeffs) {
    return(coeffs.a*x + coeffs.y0);
}

// Fit y = a*(x - x0)^b + y0
Coeffs Extrapolate::interpPowerLaw() {
    // Initialize parameters a, b, and y0
    Coeffs sol;
    double y0,y1,y2;
    double x0,x1,x2;

    sol.y0 = y0 = ys[0];
    y1 = ys[1];
    y2 = ys[2];
    sol.x0 = x0 = xs[0];
    x1 = xs[1];
    x2 = xs[2];

    // Calculate a and b coefficients, different cases to handle differing directions on the curve
    // Needs symmetry transformations to handle specific cases
    // Standard increasing x and y case (positive power)
    if ( (x2 - x0) > 0 && (x1 - x0) > 0 && (y2 - y0) > 0 && (y1 - y0) > 0 ){
        sol.b = (std::log(y2 - y0) - std::log(y1 - y0)) / (std::log(x2 - x0) - std::log(x1 - x0));
        sol.a = std::exp( std::log(y2 - y0) - sol.b*std::log(x2 - x0) );
    // Decreasing x but increasing y (negative power but backtracking)
    } else if ( (x2 - x0) < 0 && (x1 - x0) < 0 && (y2 - y0) > 0 && (y1 - y0) > 0){
        sol.b = (std::log(y2 - y0) - std::log(y1 - y0)) / (std::log(-(x2 - x0)) - std::log(-(x1 - x0)));
        sol.a = std::exp( std::log(y2 - y0) - sol.b*std::log(-(x2 - x0)) );
    // Decreasing x and decreasing y (positive power but backtrackinig)
    } else if ( (x2 - x0) < 0 && (x1 - x0) < 0 && (y2 - y0) < 0 && (y1 - y0) < 0 ){
        sol.b = (std::log( -(y2 - y0) ) - std::log( -(y1 - y0) )) / (std::log(-(x2 - x0)) - std::log(-(x1 - x0)));
        sol.a = -std::exp( std::log( -(y2 - y0) ) - sol.b*std::log(-(x2 - x0)) );
    // Increasing x but decreasing y (negative power)
    } else if ( (x2 - x0) > 0 && (x1 - x0) > 0 && (y2 - y0) < 0 && (y1 - y0) < 0 ){
        sol.b = (std::log( -(y2 - y0) ) - std::log( -(y1 - y0) )) / (std::log(x2 - x0) - std::log(x1 - x0));
        sol.a = -std::exp( std::log( -(y2 - y0) ) - sol.b*std::log(x2 - x0) );
    // Constant curve y = y0
    } else{
        sol.b = 0;
        sol.a = 0;
    }
    return sol;
}

// Fit y = a(x - x0) + y0
Coeffs Extrapolate::interpLinearLaw() {
    // Initialize parameters a, b, and c
    Coeffs sol;
    double y0,y1,y2;
    double x0,x1,x2;
    y1 = ys[1];
    y2 = ys[2];

    x1 = xs[1];
    x2 = xs[2];

    double ydiff = y2 - y1;
    double xdiff = x2 - x1;
    sol.a = ydiff / xdiff;
    sol.y0 = y2 - sol.a*(x2);

    return sol;
}

// Function to extrapolate based on the power-law fit
void Extrapolate::extrapolatePowerLaw(double x) {
    // Fit the model
    Coeffs coeffs = interpPowerLaw();

    // Predict/extrapolate the value at x
    solution = secondDerivPower(x,coeffs);
}

// Function to extrapolate based on the linear fit
void Extrapolate::extrapolateLinearLaw(double x) {
    // Fit the model
    Coeffs coeffs = interpLinearLaw();

    // Predict/extrapolate the value at x
    solution = linearModel(x,coeffs);
}

//! STILL TO DO
// Function to extrapolate based on the log fit

// constructor
Extrapolate::Extrapolate(const std::vector<double>& X, const std::vector<double>& Y) : xs(X), ys(Y) {}

