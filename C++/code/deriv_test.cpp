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
    int x_num = 101;
    int y_num = 126;
    // double x_max = -0.001;
    // double x_min = -5;
    // double y_max = 8;
    // double y_min = 1;
    double x_max = 5;
    double x_min = -5;
    double y_max = 10;
    double y_min = 5;
    std::vector<double> x,y;

    // fill grid knot arrays i.e. linspace
    // x = logspace(x_min,x_max,x_num);
    // y = logspace(y_min,y_max,y_num);
    x = linspace(x_min,x_max,x_num);
    y = linspace(y_min,y_max,y_num);
    

    // scalar field values
    // std::vector<std::vector<double>> z_gauss = gauss(x,y);
    std::vector<std::vector<double>> z_trig = trig(x,y);
    // std::vector<std::vector<double>> z_pdf = pdfLike(x,y);

    // analytic derivatives
    // Derivatives ds_gauss = deriv_gauss(x,y);
    Derivatives ds_trig = deriv_trig(x,y);
    outputTrigDerivatives(x,y,"an_ds_trig.csv");
    // Derivatives ds_pdf = deriv_pdfLike(x,y);
    // outputPdfDerivatives(x,y,"an_ds_pdf.csv");

    
    // FiniteDerivatives numds_gauss_temp(x,y,z_gauss);
    FiniteDerivatives numds_trig_temp(x,y,z_trig);
    // FiniteDerivatives numds_pdf_temp(x,y,z_pdf);

    // output derivatives to csv
    // numds_gauss_temp.outputDerivs("gauss.csv");
    numds_trig_temp.outputDerivs("num_ds_trig.csv");
    // numds_pdf_temp.outputDerivs("num_ds_pdf.csv");

    std::cout<<"DONE"<<std::endl;
    
    return 0;
    
}