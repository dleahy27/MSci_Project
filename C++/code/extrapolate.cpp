#include "../headers/extrapolate.h"

inline double square(double x){return x*x;}

// Model function y = a * x^b + c
double Extrapolate::powerLawModel(double x, const Coeffs& coeffs) {
    if ( coeffs.a*(x - coeffs.x0)<0 ){
        return (coeffs.y0 -std::pow(std::abs(coeffs.a*(x - coeffs.x0)),coeffs.b));
    }else{
        return (coeffs.y0 + std::pow((coeffs.a*(x - coeffs.x0)),coeffs.b));
    }
    
}

double Extrapolate::linearModel(double x, const Coeffs& coeffs) {
    return(coeffs.a*(x-coeffs.x0)+coeffs.y0);
}

// Fit y = a * x^b + c
Coeffs Extrapolate::interpPowerLaw() {
    // Initialize parameters a, b, and c
    Coeffs sol;
    double y0,y1,y2;
    double x0,x1,x2;
    sol.y0 = y0 = ys[0];
    y1 = ys[1];
    y2 = ys[2];
    sol.x0 = x0 = xs[0];
    x1 = xs[1];
    x2 = xs[2];

    // Step 1: Solve for b using the ratio of adjusted y values
    double log_ratio = std::log(y2 / y1);
    double x_diff_log = std::log((x2 - x0) / (x1 - x0));
    sol.b = log_ratio / x_diff_log;

    // Step 2: Solve for a using y1 and the solved b
    if ( x1-x0<0 ){
        sol.a = (y1-y0) / -std::pow(std::abs(x1 - x0), sol.b); // can also use y2 (possibly an error check in the future?)
    }else{
        sol.a = (y1-y0) / std::pow((x1 - x0), sol.b);
    }
    return sol;
}

// Gradient descent to fit y = a * x^b + c
Coeffs Extrapolate::interpLinearLaw() {
    // Initialize parameters a, b, and c
    Coeffs sol;
    double y0,y1,y2;
    double x0,x1,x2;
    sol.y0 = y0 = ys[0];
    y1 = ys[1];
    y2 = ys[2];
    sol.x0 = x0 = xs[0];
    x1 = xs[1];
    x2 = xs[2];

    double ydiff = y2 - y1;
    double xdiff = (x2 - x1);
    sol.a = ydiff / xdiff;

    return sol;
}

// Function to extrapolate based on the power-law fit
double Extrapolate::extrapolatePowerLaw(double x) {
    // Fit the model
    Coeffs coeffs = interpPowerLaw();

    // Predict/extrapolate the value at x
    return powerLawModel(x,coeffs);
}

// Function to extrapolate based on the power-law fit
double Extrapolate::extrapolateLinearLaw(double x) {
    // Fit the model
    Coeffs coeffs = interpLinearLaw();

    // Predict/extrapolate the value at x
    return linearModel(x,coeffs);
}

// constructor
Extrapolate::Extrapolate(double x, const std::vector<double>& X, const std::vector<double>& Y) : xs(X), ys(Y) {
    solution = extrapolatePowerLaw(x);
    //solution = extrapolateLinearLaw(x);
}

