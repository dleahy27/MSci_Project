#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include "cubicspline.h"


std::vector<double> gauss_trig(const std::vector<double>& x){
     int n = x.size();
    
    std::vector<double> y(n);
    
    for (  int i = 0; i<=n; i++ ){
        y[i] = (std::cos(3*x[i]) + 1.5) * std::exp(-(std::pow((x[i]-3),2)/6));
    }
    
    return y;
}

double gauss_trig(const double x){
    return (std::cos(3*x) + 1.5) * std::exp(-(std::pow((x-3),2)/6));
}

int main(int argc, char* argv[]){
    int num = 100;
    double max = 15;
    double min = -15;
    std::vector<double> x;

    for (  int i = 0; i<num; i++){
        double val = min + i*(max - min)/num; 
        x.push_back(val);
    }

    std::vector<double> y = gauss_trig(x);

    CubicSpline spline(x,y,0,0);

    std::vector<std::vector<double>> sol = spline.GetS();
    
    // std::cout << "The spline coefficients of each interval are " << std::endl;
    // for ( const auto& row : sol ){
    //     for ( const auto& val : row ){
    //         std::cout << val << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // list of values test
    int num_t = 10;
    double max_t = 12.456;
    double min_t = -12.456;
    std::vector<double> x_ts;
    for (  int i = 0; i<num_t; i++){
        double val = min_t + i*(max_t - min_t)/num_t; 
        x_ts.push_back(val);
    }
    std::vector<double> y_ts = gauss_trig(x_ts);
    std::vector<double> y_sols = spline.evaluateSpline(x_ts, sol, x, true);

    //single value test
    double x_t = 0.69420;
    double y_t = gauss_trig(x_t);
    double y_sol = spline.evaluateSpline(x_t, sol, x, true);
    

    std::cout << "Function value at the point "<< x_t <<" is:" << std::endl;
    std::cout<<" x  |  y  |  y_s " <<  std::endl<< "------------------"<<std::endl;
    std::cout<<x_t<<" "<<y_t << " " << y_sol << std::endl;
    
    std::cout << "Function values at "<< num_t << " equally spaced points between " << max_t << " and " << min_t << " are:" << std::endl;
    std::cout<<" x  |  y  |  y_s " <<  std::endl<< "------------------"<<std::endl;
    for (  int i = 0; i<num_t; i++){
        std::cout<<x_ts[i]<<" "<<y_ts[i] << " " << y_sols[i] << std::endl;
    }

    return 0;
}