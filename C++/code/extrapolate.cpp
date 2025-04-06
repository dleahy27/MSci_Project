#include "../headers/extrapolate.h"

// fast square function
inline double square(double t){return t*t;}
inline const double epsilon = std::numeric_limits<double>::epsilon();

// Model function y = a(t - t0)^b + xf0
double Extrapolate::powerLawModel(double t, const Coeffs& coeffs) {
    // Covers some backwards cases by reflecting around t-t0
    if ( (t - coeffs.t0)<0 ){
        return (coeffs.xf0 + coeffs.a*std::pow((-t + coeffs.t0),coeffs.b));
    }else{
        return (coeffs.xf0 + coeffs.a*std::pow(((t - coeffs.t0)),coeffs.b));
    }
}

// Analytical second derivative of power function
double Extrapolate::secondDerivPower(double t, const Coeffs& coeffs) {
    // Covers some backwards cases by reflecting around t-t0
    if (std::abs(coeffs.b-1)<=100*epsilon){
        return 0;
    }else if ( t - coeffs.t0<0 ){
        return coeffs.a*coeffs.b *(coeffs.b - 1)*std::pow(-t + coeffs.t0, coeffs.b - 2);
    }else if ( t - coeffs.t0>0 ){
        return  coeffs.a*coeffs.b *(coeffs.b - 1)*std::pow(t - coeffs.t0, coeffs.b - 2);
    }else {
        return 0;
    }
}

// Model function y = a (10^(t - t0)^b) + xf0
double Extrapolate::tenPowerLawModel(double t, const Coeffs& coeffs) {
    return coeffs.a*std::pow(10,t*coeffs.b) - coeffs.xf0;
}

// Analytical second derivative of power function
double Extrapolate::secondDerivTenPower(double t, const Coeffs& coeffs) {
    if (std::abs(coeffs.b)<=100*epsilon){
        return 0;
    }else{return square(std::log(10)*coeffs.b)*coeffs.a*std::pow(10, t*coeffs.b);}
}

// Linear model: y = a(t-t0) + xf0
double Extrapolate::linearModel(double t, const Coeffs& coeffs) {
    return(coeffs.a*t + coeffs.xf0);
}

// Fit y = a*(t - t0)^b + xf0
Coeffs Extrapolate::interpPowerLaw() {
    // Initialize parameters a, b, and xf0
    Coeffs sol;
    double xf0,xf1,xf2;
    double t0,t1,t2;

    sol.xf0 = xf0 = xfs[0];
    xf1 = xfs[1];
    xf2 = xfs[2];
    sol.t0 = t0 = ts[0];
    t1 = ts[1];
    t2 = ts[2];

    // Calculate a and b coefficients, different cases to handle differing directions on the curve
    // Needs symmetry transformations to handle specific cases
    // Standard increasing t and y case (positive power)
    if ( (t2 - t0) > 0 && (t1 - t0) > 0 && (xf2 - xf0) > 0 && (xf1 - xf0) > 0 ){
        sol.b = (std::log(xf2 - xf0) - std::log(xf1 - xf0)) / (std::log(t2 - t0) - std::log(t1 - t0));
        sol.a = std::exp( std::log(xf2 - xf0) - sol.b*std::log(t2 - t0) );
        return sol;
        // Decreasing t but increasing y (negative power but backtracking)
    } else if ( (t2 - t0) < 0 && (t1 - t0) < 0 && (xf2 - xf0) > 0 && (xf1 - xf0) > 0){
        sol.b = (std::log(xf2 - xf0) - std::log(xf1 - xf0)) / (std::log(-(t2 - t0)) - std::log(-(t1 - t0)));
        sol.a = std::exp( std::log(xf2 - xf0) - sol.b*std::log(-(t2 - t0)) );
        return sol;
        // Decreasing t and decreasing y (positive power but backtrackinig)
    } else if ( (t2 - t0) < 0 && (t1 - t0) < 0 && (xf2 - xf0) < 0 && (xf1 - xf0) < 0 ){
        sol.b = (std::log( -(xf2 - xf0) ) - std::log( -(xf1 - xf0) )) / (std::log(-(t2 - t0)) - std::log(-(t1 - t0)));
        sol.a = -std::exp( std::log( -(xf2 - xf0) ) - sol.b*std::log(-(t2 - t0)) );
        return sol;
        // Increasing t but decreasing y (negative power)
    } else if ( (t2 - t0) > 0 && (t1 - t0) > 0 && (xf2 - xf0) < 0 && (xf1 - xf0) < 0 ){
        sol.b = (std::log( -(xf2 - xf0) ) - std::log( -(xf1 - xf0) )) / (std::log(t2 - t0) - std::log(t1 - t0));
        sol.a = -std::exp( std::log( -(xf2 - xf0) ) - sol.b*std::log(t2 - t0) );
        return sol;
        // Constant curve y = xf0
    } else{
        sol.b = 0;
        sol.a = 0;
        return sol;
    }
}

// Fit y = a(t - t0) + xf0
Coeffs Extrapolate::interpLinearLaw() {
    // Initialize parameters a, b, and c
    Coeffs sol;
    double xf0,xf1,xf2;
    double t0,t1,t2;
    xf1 = xfs[1];
    xf2 = xfs[2];

    t1 = ts[1];
    t2 = ts[2];

    double ydiff = xf2 - xf1;
    double xdiff = t2 - t1;
    sol.a = ydiff / xdiff;
    sol.xf0 = xf2 - sol.a*(t2);

    return sol;
}

// Fit y = a* (10^(t - t0))^b + xf0
Coeffs Extrapolate::interpTenPowerLaw() {
    // Initialize parameters a, b, and xf0
    Coeffs sol;
    double xf1,xf2;
    double t1,t2;

    t1 = ts[0];
    t2 = ts[1];

    if (xfs.size() == 3){
        sol.xf0 = xfs[0];

        if (xfs[2] - xfs[0] > 0 && xfs[1] - xfs[0] > 0){
            xf1 = std::log10(xfs[1] - xfs[0]);
            xf2 = std::log10(xfs[2] - xfs[0]);
            sol.b = (xf2 - xf1) / (t2 - t1);
            sol.a = std::pow(10,(xf2 - sol.b*t2));

            return sol; 
        } else if (xfs[2] - xfs[0] < 0 && xfs[1] - xfs[0] < 0){
            xf1 = std::log10(xfs[0] - xfs[1]);
            xf2 = std::log10(xfs[0] - xfs[2]);
            sol.b = (xf2 - xf1) / (t2 - t1);
            sol.a = std::pow(10,(xf2 - sol.b*t2));

            return sol;
        } else{
            sol.a = 0;
            sol.b = 0;

            return sol;
        }
    } else if (xfs.size() == 2){
        sol.xf0 = 0;

        if (xfs[1]> 0 && xfs[0]> 0){
            xf1 = std::log10(xfs[0]);
            xf2 = std::log10(xfs[1]);
            sol.b = (xf2 - xf1) / (t2 - t1);
            sol.a = std::pow(10,(xf2 - sol.b*t2));
            return sol; 
        } else if (xfs[2] - xfs[0] < 0 && xfs[1] - xfs[0] < 0){
            xf1 = std::log10(-xfs[0]);
            xf2 = std::log10(-xfs[1]);
            sol.b = (xf2 - xf1) / (t2 - t1);
            sol.a = std::pow(10,(xf2 - sol.b*t2));
            return sol;
        } else{
            sol.a = 0;
            sol.b = 0;

            return sol;
        }
    }else {
        // RETURN ERROR
        return sol;
    }
}

// Function to extrapolate based on the power-law fit
void Extrapolate::extrapolatePowerLaw(double t) {
    // Fit the model
    Coeffs coeffs = interpPowerLaw();

    // Predict/extrapolate the value at t
    solution = secondDerivPower(t,coeffs);
}

// Function to extrapolate based on the 10^power-law fit
void Extrapolate::extrapolateTenPowerLaw(double t) {
    // Fit the model
    Coeffs coeffs = interpTenPowerLaw();

    // Predict/extrapolate the value at t
    solution = secondDerivTenPower(t,coeffs);
}

// Function to extrapolate based on the linear fit
void Extrapolate::extrapolateLinearLaw(double t) {
    // Fit the model
    Coeffs coeffs = interpLinearLaw();

    // Predict/extrapolate the value at t
    solution = linearModel(t,coeffs);
}

// constructor
Extrapolate::Extrapolate(const std::vector<double>& X, const std::vector<double>& Y) : ts(X), xfs(Y) {}

double Extrapolate::powerLawxfs(double t){
    // Fit the model
    Coeffs coeffs = interpPowerLaw();

    // Predict/extrapolate the value at t
    return(powerLawModel(t,coeffs));
}

double Extrapolate::tenPowerLawxfs(double t){
    // Fit the model
    Coeffs coeffs = interpTenPowerLaw();

    // Predict/extrapolate the value at t
    return(tenPowerLawModel(t,coeffs));

}
