#include "../headers/boundaryderivatives.h"

// fast square function
inline double square(double x){return x*x;}

// fast integer power function
inline constexpr double power(double base, int exponent) {
    if (exponent == 0) {
        return 1;
    } else if (exponent > 0) {
        return base * power(base, exponent - 1);
    } else {
        return std::pow(base,exponent);
    }
}

// constructor
BoundaryDerivatives::BoundaryDerivatives(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<std::vector<double>>& Z) : x(X), y(Y), z(Z), m(x.size()), n(y.size()), d2_left(n), d2_right(n), d2_bottom(m), d2_top(m), d4y2x2(4) {
    // calculate boundary derivatives
    partialDerivative2Y();
    partialDerivative2X();
    partialDerivative2X2Y();
    // set them to struct
    SetBoundaryDerivatives();
}

// Second y partial
void BoundaryDerivatives::partialDerivative2Y() {
    // Extrapolate to boundary 
    std::vector<double> tempVec1 = {y[3],y[2],y[1]};
    std::vector<double> tempVec3 = {y[n-4],y[n-3],y[n-2]};
    for (int i = 0; i < m; i++) {
        std::vector<double> tempVec2 = {z[3][i],z[2][i],z[1][i]};
        std::vector<double> tempVec4 = {z[n-4][i],z[n-3][i],z[n-2][i]};

        Extrapolate fit1(tempVec1, tempVec2);
        fit1.extrapolatePowerLaw(y[0]);
        Extrapolate fit2(tempVec3, tempVec4);
        fit2.extrapolatePowerLaw(y[n-1]);

        //std::cout<<fit1.solution<<std::endl;
        d2_bottom[i] = fit1.solution;
        d2_top[i] = fit2.solution;
    }
}

// Second x partial
void BoundaryDerivatives::partialDerivative2X() {
    // Extrapolate to boundary 
    std::vector<double> tempVec1 = {x[3],x[2],x[1]};
    std::vector<double> tempVec3 = {x[m-4],x[m-3],x[m-2]};
    for (int j = 0; j < n; j++) {
        std::vector<double> tempVec2 = {z[j][3],z[j][2],z[j][1]};
        std::vector<double> tempVec4 = {z[j][m-4],z[j][m-3],z[j][m-2]};
        
        Extrapolate fit1(tempVec1, tempVec2);
        fit1.extrapolatePowerLaw(x[0]);
        Extrapolate fit2(tempVec3, tempVec4);
        fit2.extrapolatePowerLaw(x[m-1]);

        //std::cout<<fit1.solution<<std::endl;
        d2_left[j] = fit1.solution;
        d2_right[j] = fit2.solution;
    }
}

void BoundaryDerivatives::partialDerivative2X2Y() {
    // Assume x and y are independant, d4 is then product of d2s
    // corner derivatives go counter-clockwise: [0,0] [1,0] [1,1] [0,1]
    // Put in internal variable
    d4y2x2[0] = d2_left[0]*d2_bottom[0];
    d4y2x2[1] = d2_right[0]*d2_bottom[m-1];
    d4y2x2[2] = d2_right[n-1]*d2_top[m-1];
    d4y2x2[3] = d2_left[n-1]*d2_top[0];
}


// Sets the bounadry derivatives to the struct
void BoundaryDerivatives::SetBoundaryDerivatives(){
    boundary_derivatives.d2x2s_left = d2_left;
    boundary_derivatives.d2y2s_top = d2_top;
    boundary_derivatives.d2x2s_right = d2_right;
    boundary_derivatives.d2y2s_bottom = d2_bottom;
    boundary_derivatives.d4x2y2s_corners = d4y2x2;
    
}

// Testing function that writes out the derivative values
void BoundaryDerivatives::outputDerivs(std::string filename){
    outputXDerivs(filename);
    outputYDerivs(filename); 
    outputCornerDerivs(filename);  
}

void BoundaryDerivatives::outputYDerivs(std::string filename){
    std::ofstream myfile;
    myfile.open("../outputs/y_"+filename);
    myfile << "x,d2y_top,d2y_bottom"<<std::endl;
    for(int i = 0; i<m; i++){
        myfile<<x[i]<<","<<d2_top[i]<<","<<d2_bottom[i]<<std::endl;
    }
    myfile.close();
}

void BoundaryDerivatives::outputXDerivs(std::string filename){
    std::ofstream myfile;
    myfile.open("../outputs/x_"+filename);
    myfile << "y,d2x_left,d2x_right"<<std::endl;
    for(int i = 0; i<n; i++){
        myfile<<y[i]<<","<<d2_left[i]<<","<<d2_right[i]<<std::endl;
    }
    myfile.close();
}

void BoundaryDerivatives::outputCornerDerivs(std::string filename){
    std::ofstream myfile;
    myfile.open("../outputs/corner_"+filename);
    myfile << "x,y,d4d2xd2y"<<std::endl;
    myfile<<x[0]<<","<<y[0]<<","<<d4y2x2[0]<<std::endl;
    myfile<<x[m]<<","<<y[0]<<","<<d4y2x2[1]<<std::endl;
    myfile<<x[m]<<","<<y[n]<<","<<d4y2x2[2]<<std::endl;
    myfile<<x[0]<<","<<y[n]<<","<<d4y2x2[3]<<std::endl;
    myfile.close();
}