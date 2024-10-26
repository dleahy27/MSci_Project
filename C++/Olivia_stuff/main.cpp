#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include "bicubicspline.h"


std::vector<std::vector<double>> gauss(std::vector<double> x, std::vector<double> y){
    int n,m;

    m = x.size();
    n = y.size();
    
std::vector<std::vector<double>> z(m, std::vector<double> (n,0));
    
    for (  int i=0; i<m; i++ ){
        for (  int j=0; j<n; j++  ){
            z[i][j] = (std::cos(3*y[j]) * std::sin(3*x[i]) + 1.5) * std::exp(-(std::pow((x[i]-3),2) + std::pow((y[j]-3),2)) / 6);
        }
    }
    
    return z;
}

double gauss(double x, double y){return (std::cos(3*y) * std::sin(3*x) + 1.5) * std::exp(-(std::pow((x-3),2) + std::pow((y-3),2)) / 6);}

int main(int argc, char* argv[]){
    // definitions
    int num = 100;
    double max = 5;
    double min = -5;
    std::vector<double> x,y;

    // fill grid knot arrays
    for (  int i = 0; i<num; i++){
        double val = min + i*(max - min)/num; 
        x.push_back(val);
        y.push_back(val);
    }
    
   
    // scalar field values
    std::vector<std::vector<double>> z = gauss(x,y);

    for ( const auto& row : z ){
        for (const auto& val : row ){
            std::cout<< val << " ";
        }
        std::cout<<std::endl;
    }

    // initialisation of class test
    BicubicSpline spline(x,y,z);

    

    // single value test
    double x_t = 1.69;
    double y_t = 4.201;

    double z_t = gauss(x_t,y_t);
    double z_spline = spline.evaluateSpline(x_t,y_t);

    std::cout << "Function value at the point ("<< x_t << y_t <<") is:" << std::endl;
    std::cout<<" x  |  y  |  z  |  z_s " <<  std::endl<< "------------------------"<<std::endl;
    std::cout<<x_t<<" "<<y_t << " " << z_t << " " << z_spline << std::endl;
    

    // array of values test
    int y_n = 3;
    int x_n = 3;
    double max_x = 3.6598;
    double min_x = -3.6598;
    double max_y = 4.887;
    double min_y = -4.887;

    std::vector<double> xs_t, ys_t;
    std::vector<std::vector<double>> zs_t,zs_spline;

    for (  int i = 0; i<x_n; i++){xs_t.push_back(min_x + i*(max_x - min_x)/x_n);}
    for (  int i = 0; i<y_n; i++){ys_t.push_back(min_y + i*(max_y - min_y)/y_n);}

    zs_t = gauss(xs_t,ys_t);
    zs_spline = spline.evaluateSpline(xs_t,ys_t);

    std::cout << "Function values at ["<< x_n << y_n << "] equally spaced points on the grid [" << max_x << "," << min_x << "]x["<< max_y << "," << min_y << "] are:" << std::endl;
    std::cout<<" x  |  y  |  y_s " <<  std::endl<< "------------------"<<std::endl;
    for (  int i = 0; i<x_n; i++){
        for (  int j = 0; j<y_n; j++){
            std::cout<<xs_t[i]<<" "<<ys_t[j] << " " << zs_t[i][j] << " " << zs_spline[i][j] << std::endl;
        }
    }

    return 0;
}