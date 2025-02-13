#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include "../headers/bicubicspline.h"
#include "../headers/boundaryderivatives.h"
#include "../headers/funcs.h"
#include "../headers/sample.h"
#include "../headers/derivatives.h"

using namespace funcs;
using namespace sample;

void deleteDirContent(const std::filesystem::path& dir_path) {
    for (auto& path: std::filesystem::directory_iterator(dir_path)) {
        std::filesystem::remove_all(path);
    }
}

void pdfTest(int x_num, int y_num, int n){
    double x_max = -0.001;
    double x_min = -5;
    double y_max = 8;
    double y_min = 1;
    std::vector<double> x,y;

    // fill grid knot arrays i.e. linspace
    x = linspace(x_min,x_max,x_num);
    y = linspace(y_min,y_max,y_num);
    
    // scalar field values
    std::vector<std::vector<double>> z = pdfLike(x,y);

    // numerical derivatives
    Derivatives ds_num;
    BoundaryDerivatives ds_num_temp(x,y,z);
    ds_num = ds_num_temp.boundary_derivatives;

    // initialisation of spline);
    BicubicSpline spline_ds_num(x,y,z,ds_num);


    // array of values test
    // initialisations
    int y_n = 101;
    int x_n = 101;
    double max_x = -0.0011;
    double min_x = -4.99;
    double max_y = 7.99;
    double min_y = 1.01;

    std::vector<double> xs_t, ys_t;
    std::vector<std::vector<double>> zs_t, zs_spline_ds_num;

    // set xy vectors i.e. linspace
    xs_t = linspace(min_x, max_x, x_n);
    ys_t = linspace(min_y, max_y, y_n);
    
    // evaluate splines + function
    zs_t = pdfLike(xs_t,ys_t);
    zs_spline_ds_num = spline_ds_num.evaluateSpline(xs_t,ys_t);
    // remove old file then write the data to file
    std::filesystem::remove("../outputs/grid_x_test/pdf_xgrid_"+ std::to_string(n) +".csv");
    std::ofstream myfile;
    myfile.open("../outputs/grid_x_test/pdf_xgrid_"+ std::to_string(n) +".csv");
    myfile << "nx,ny,x,y,z,z_spline\n";
    for (  int i = 0; i<y_n; i++){
        for (  int j = 0; j<x_n; j++){
            myfile<<x_num<<","<<y_num<<","<<xs_t[j]<<","<<ys_t[i] << "," << zs_t[i][j] << "," <<zs_spline_ds_num[i][j] << std::endl;
        }
    }
    myfile.close();
}
int main(int argc, char* argv[]){
    // remove old files
    deleteDirContent("../outputs/grid_x_test/");
    // Initialisations
    std::vector<int> x_num = linspace(21,71,51);
    std::vector<int> y_num = {31};

    int n = 0;
    for ( auto const& x : x_num ){
        for ( auto const& y : y_num ){
            pdfTest(x,y,n);
            n++;
        }
    }
    
    std::cout<<"DONE"<<std::endl;

    return 0;
} 