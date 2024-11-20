#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "../headers/bicubicspline.h"
#include "../headers/finitederivatives.h"
#include "../headers/funcs.h"
#include "../headers/derivatives.h"

using namespace funcs;

int main(int argc, char* argv[]){
    // Initialisations
    int x_num = 51;
    int y_num = 51;
    double x_max = 5;
    double x_min = -5;
    double y_max = 5;
    double y_min = -5;
    std::vector<double> x,y;

    
    // fill grid knot arrays i.e. linspace
    x = initialisation(x_min,x_max,x_num);
    y = initialisation(y_min,y_max,y_num);
    
    // scalar field values
    std::vector<std::vector<double>> z = gauss(x,y);

    // analytic derivatives
    Derivatives ds = deriv_gauss(x,y);

    // numerical derivatives
    Derivatives ds_num;
    FiniteDerivatives ds_num_temp(x,y,z);
    ds_num = ds_num_temp.boundary_derivatives;
    // initialisation of class test both with derivaties and without
    BicubicSpline spline(x,y,z);
    BicubicSpline spline_ds(x,y,z,ds);
    BicubicSpline spline_ds_num(x,y,z,ds_num);

    // array of values test
    // initialisations
    int y_n = 27;
    int x_n = 27;
    double max_x = 4.99;
    double min_x = -4.99;
    double max_y = 4.99;
    double min_y = -4.99;

    std::vector<double> xs_t, ys_t;
    std::vector<std::vector<double>> zs_t, zs_spline, zs_spline_ds_num, zs_spline_ds;

    // set xy vectors i.e. linspace
    xs_t = initialisation(min_x, max_x, x_n);
    ys_t = initialisation(min_y, max_y, y_n);
    
    // evaluate splines + function
    zs_t = gauss(xs_t,ys_t);
    zs_spline_ds = spline_ds.evaluateSpline(xs_t,ys_t);
    zs_spline_ds_num = spline_ds_num.evaluateSpline(xs_t,ys_t);
    zs_spline = spline.evaluateSpline(xs_t,ys_t);

    // write the data to file
    std::ofstream myfile;
    myfile.open("../outputs/gauss_bicubic.csv");
    myfile << "x,y,z,z_spline_ds,z_spline_ds_num,z_spline\n";
    for (  int i = 0; i<y_n; i++){
        for (  int j = 0; j<x_n; j++){
            myfile<<xs_t[j]<<","<<ys_t[i] << "," << zs_t[i][j] <<"," <<zs_spline_ds[i][j] << "," <<zs_spline_ds_num[i][j] << "," << zs_spline[i][j] << std::endl;
        }
    }
    myfile.close();
    std::cout<<"DONE"<<std::endl;
    return 0;
}