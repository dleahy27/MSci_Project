#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include "../headers/bicubicspline.h"
#include "../headers/finitederivatives.h"
#include "../headers/funcs.h"
#include "../headers/derivatives.h"
#include "../headers/sample.h"

using namespace funcs;
using namespace sample;

double splineError(const std::vector<std::vector<double>>& z1, const std::vector<std::vector<double>>& z2){
    double err;
    for (int i = 0; i<z1.size(); i++){
        for (int j = 0; j<z1[0].size(); j++){
            err += std::abs(z1[i][j]-z2[i][j]);
        }
    }
    return err;
}

std::vector<std::vector<double>> timeTest(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& z, const std::vector<double>& xs_t, const std::vector<double>& ys_t, const Derivatives& ds_num){
    std::cout<<"13"<<std::endl;
    // initialisation of class 
    BicubicSpline spline_ds_num(x,y,z,ds_num);
    std::cout<<"14"<<std::endl;
    // evaluate splines at array of values
    return{spline_ds_num.evaluateSpline(xs_t,ys_t)};

}

int main(int argc, char* argv[]){
    // Initialisations
    std::vector<int> ms,ns;
    double x_max = 5;
    double x_min = -5;
    double y_max = 5;
    double y_min = -5;

    int m_min = 50;
    int n_min = 50;
    int m_max = 151;
    int n_max = 151;

    ms = linspace(m_min, m_max, 101);
    ns = linspace(n_min, n_max, 101);
    // array of values test
    // initialisations
    int y_n = 27;
    int x_n = 27;
    double max_x = 4.99;
    double min_x = -4.99;
    double max_y = 4.99;
    double min_y = -4.99;
    std::vector<double> xs_t, ys_t;

    // set xy vectors i.e. linspace
    xs_t = linspace(min_x, max_x, x_n);
    ys_t = linspace(min_y, max_y, y_n); 
    
    std::ofstream myfile;
    myfile.open("../outputs/accuracy_test.csv");
    myfile << "m,n,time,error\n";

    for (int i = 0; i<10; i++){
        //for (int j = 0; j<ns.size(); j++){
        myfile<<ms[i]<<","<<ns[i]<<",";

        std::vector<double> x = linspace(x_min,x_max,ms[i]);
        std::vector<double> y = linspace(y_min,y_max,ns[i]);

        std::vector<std::vector<double>>z = gauss(x,y);

        FiniteDerivatives ds_num_temp(x,y,z);
        Derivatives ds_num = ds_num_temp.boundary_derivatives;

        const auto start = std::chrono::steady_clock::now();
        std::vector<std::vector<double>> z_num = timeTest(x,y,z,xs_t,ys_t,ds_num);
        const auto end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = end - start;
        myfile<<elapsed.count()<<",";

        double error = splineError(z,z_num);
        myfile<<error<<std::endl;
    }

    myfile.close();
    
    return 0;
}