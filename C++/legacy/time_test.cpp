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

inline void timeTest(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::vector<double>>& z, const std::vector<double>& xs_t, const std::vector<double>& ys_t, const Derivatives& ds_num){
    

}

int main(int argc, char* argv[]){
    int N = 1000;

    // Initialisations
    int x_num = 51;
    int y_num = 51;
    double x_max = 5;
    double x_min = -5;
    double y_max = 5;
    double y_min = -5;
    std::vector<double> x,y;

    // fill grid knot arrays i.e. linspace
    x = linspace(x_min,x_max,x_num);
    y = linspace(y_min,y_max,y_num);

    // array of values test
    // initialisations
    int y_n = 100;
    int x_n = 100;

    std::vector<double> xs_t, ys_t;

    // set xy vectors i.e. linspace
    xs_t = random_uniform(x_min, x_max, x_n);
    ys_t = random_uniform(y_min, y_max, y_n);

    // scalar field values
    std::vector<std::vector<double>> z = trig(x,y);

    // numerical derivatives
    Derivatives ds_num;
    FiniteDerivatives ds_num_temp(x,y,z);
    ds_num = ds_num_temp.boundary_derivatives;

    std::ofstream myfile;
    myfile.open("../outputs/time_test_1000.csv");
    myfile << "t_init,t_eval\n";

    for (int i = 0; i<N; i++){
        auto start = std::chrono::steady_clock::now();
        // initialisation of class 
        BicubicSpline spline_ds_num(x,y,z,ds_num);
        auto end = std::chrono::steady_clock::now();                  \
        const std::chrono::duration<double> elapsed1 = end - start;
        // evaluate splines at array of values
        start = std::chrono::steady_clock::now();
        spline_ds_num.evaluateSpline(xs_t,ys_t);
        end = std::chrono::steady_clock::now();                  \
        const std::chrono::duration<double> elapsed2 = end - start; 
        myfile<<elapsed1.count()<<","<<elapsed2.count()<<std::endl;
    }

    myfile.close();
    std::cout<<"DONE"<<std::endl;
    return 0;
}