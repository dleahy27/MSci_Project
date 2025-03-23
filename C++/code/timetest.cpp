#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include "../headers/bicubicspline.h"
#include "../headers/boundaryderivatives.h"
#include "../headers/funcs.h"
#include "../headers/sample.h"
#include "../headers/derivatives.h"

using namespace funcs;
using namespace sample;

int main(int argc, char* argv[]){
    std::vector<double> y = {1.300000e+00,1.501480e+00,1.754940e+00,2.077770e+00,2.494520e+00,3.040400e+00,3.766790e+00,4.750000e+00,6.133070e+00,8.090760e+00,1.092470e+01,1.512840e+01,2.153060e+01,3.156460e+01,4.778650e+01,7.491100e+01,1.219530e+02,2.068370e+02,3.667270e+02,6.822740e+02,1.337310e+03,2.773760e+03,6.116880e+03,1.441630e+04,3.651470e+04,1.000000e+05};
    std::vector<double> x = {1.000000e-08,1.214290e-08,1.474520e-08,1.790520e-08,2.174240e-08,2.640200e-08,3.206010e-08,3.893070e-08,4.727370e-08,5.740470e-08,6.970670e-08,8.464500e-08,1.027840e-07,1.248110e-07,1.515580e-07,1.840360e-07,2.234750e-07,2.713640e-07,3.295150e-07,4.001270e-07,4.858690e-07,5.899830e-07,7.164050e-07,8.699160e-07,1.056320e-06,1.282650e-06,1.557480e-06,1.892460e-06,2.297910e-06,2.790210e-06,3.387960e-06,4.113730e-06,4.995000e-06,6.065420e-06,7.365210e-06,8.943510e-06,1.086000e-05,1.318710e-05,1.601280e-05,1.944390e-05,2.361000e-05,2.866850e-05,3.481040e-05,4.226760e-05,5.132150e-05,6.231370e-05,7.565850e-05,9.185860e-05,1.115240e-04,1.353930e-04,1.643640e-04,1.995060e-04,2.421600e-04,2.939080e-04,3.566770e-04,4.327950e-04,5.250750e-04,6.369110e-04,7.726750e-04,9.367690e-04,1.135340e-03,1.375450e-03,1.665550e-03,2.015680e-03,2.437770e-03,2.945870e-03,3.556390e-03,4.288450e-03,5.164110e-03,6.208550e-03,7.450380e-03,8.921250e-03,1.065330e-02,1.268820e-02,1.506290e-02,1.781700e-02,2.098570e-02,2.461440e-02,2.873900e-02,3.338270e-02,3.857250e-02,4.433070e-02,5.066750e-02,5.758780e-02,6.509020e-02,7.316730e-02,8.180630e-02,9.099030e-02,1.006990e-01,1.109100e-01,1.215990e-01,1.327420e-01,1.443130e-01,1.562890e-01,1.686470e-01,1.813520e-01,1.944020e-01,2.077620e-01,2.214160e-01,2.353440e-01,2.495390e-01,2.639710e-01,2.786310e-01,2.935040e-01,3.085780e-01,3.238410e-01,3.392820e-01,3.548990e-01,3.706620e-01,3.865750e-01,4.026290e-01,4.188160e-01,4.351300e-01,4.515630e-01,4.681100e-01,4.847630e-01,5.015190e-01,5.183720e-01,5.353160e-01,5.523490e-01,5.694650e-01,5.866600e-01,6.039320e-01,6.212820e-01,6.386840e-01,6.561650e-01,6.737080e-01,6.913110e-01,7.089700e-01,7.266830e-01,7.444490e-01,7.622630e-01,7.801240e-01,7.980270e-01,8.159700e-01,8.339490e-01,8.519580e-01,8.699910e-01,8.880430e-01,9.060860e-01,9.240980e-01,9.420120e-01,9.596630e-01,9.764990e-01,9.901900e-01,9.965610e-01,9.985400e-01,9.993360e-01,9.997490e-01,1.000000e+00};
    
    for (uint i = 0; i<x.size(); i++){
        x[i] = std::log10(x[i]);
    }
    for (uint j = 0; j<y.size(); j++){
        y[j] = 2*std::log10(y[j]);
    }
    std::ifstream file("../test_down_xfx.dat"); // Open the file
    std::string line;
    std::vector<std::vector<double>> z(y.size(), std::vector<double>(x.size())); // 2D vector to hold the data
    std::vector<double> temp;

    // Read each line from the file
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            temp.push_back(std::stod(value));
        }

        for (uint i = 0; i<x.size(); i++){
            for (uint j = 0; j<y.size(); j++){
                z[j][i] = temp[i*y.size()+ j];
            }
        }
    }

    // write the data to file
    std::ofstream myfile;
    myfile.open("../outputs/times.csv");
    myfile << "x_n,y_n,t\n";
    // numerical derivatives
    auto start = std::chrono::steady_clock::now();
    Derivatives ds_num;
    BoundaryDerivatives ds_num_temp(x,y,z);
    ds_num = ds_num_temp.boundary_derivatives;

    BicubicSpline spline_ds_num(x,y,z);
    auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed = end - start;
    myfile <<0<<","<<0<<","<<elapsed.count()<< std::endl;
    // array of values test
    // initialisations
    int y_n = 1000;
    std::vector<int> x_n = linspace(101,10001,1000);
    double max_x = 0;
    double min_x = -8;
    double max_y = 10;
    double min_y = 1;

    std::vector<double> xs_t, ys_t;
    std::vector<std::vector<double>> zs_spline_ds_num;

    // set xy vectors i.e. linspace
    
    ys_t = linspace(min_y, max_y, y_n);
    
    for(int i = 0; i<x_n.size(); i++){
        // evaluate splines + function
        xs_t = linspace(min_x, max_x, x_n[i]);
        auto start = std::chrono::steady_clock::now();
        spline_ds_num.evaluateSpline(xs_t,ys_t);
        auto end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = end - start;
        myfile<<y_n<<","<<x_n[i]<<","<<elapsed.count()<< std::endl;
    }
    myfile.close();
    std::cout<<"DONE"<<std::endl;
    return 0;
}