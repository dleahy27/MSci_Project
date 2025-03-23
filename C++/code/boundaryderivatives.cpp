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
BoundaryDerivatives::BoundaryDerivatives(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<std::vector<double>>& Z) : x(X), y(Y), z(Z), m(x.size()), n(y.size()), d2_left(n), d2_right(n), d2_bottom(m), d2_top(m), d4y2x2(4), z_left(n), z_right(n), z_bottom(m), z_top(m) {
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
    std::vector<double> tempVec1 = {y[2],y[1],y[0]};
    std::vector<double> tempVec3 = {y[n-3],y[n-2],y[n-1]};
    for (int i = 0; i < m; i++) {
        std::vector<double> tempVec2 = {z[2][i],z[1][i],z[0][i]};
        std::vector<double> tempVec4 = {z[n-3][i],z[n-2][i],z[n-1][i]};

        Extrapolate fit1(tempVec1, tempVec2);
        fit1.extrapolatePowerLaw(y[0]);
        Extrapolate fit2(tempVec3, tempVec4);
        fit2.extrapolatePowerLaw(y[n-1]);

        d2_bottom[i] = fit1.solution;
        d2_top[i] = fit2.solution;

        z_bottom[i] = fit1.powerLawZs(y[0]);
        z_top[i] = fit2.powerLawZs(y[n-1]);
    }
}

// Second x partial
void BoundaryDerivatives::partialDerivative2X() {
    // Extrapolate to boundary 
    std::vector<double> tempVec1 = {x[1],x[0]};
    std::vector<double> tempVec3 = {x[m-2],x[m-1]};
    for (int j = 0; j < n; j++) {
        std::vector<double> tempVec2 = {z[j][1],z[j][0]}; 
        std::vector<double> tempVec4 = {z[j][m-3],z[j][m-2],z[j][m-1]};
        
        Extrapolate fit1(tempVec1, tempVec2);
        fit1.extrapolateTenPowerLaw(x[0]);
        Extrapolate fit2(tempVec3, tempVec4);
        fit2.extrapolateTenPowerLaw(x[m-1]);

        d2_left[j] = fit1.solution;
        d2_right[j] = fit2.solution;

        z_left[j] = fit1.tenPowerLawZs(x[0]);
        z_right[j] = fit2.tenPowerLawZs(x[m-1]);
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
    myfile<<x[m-1]<<","<<y[0]<<","<<d4y2x2[1]<<std::endl;
    myfile<<x[m-1]<<","<<y[n-1]<<","<<d4y2x2[2]<<std::endl;
    myfile<<x[0]<<","<<y[n-1]<<","<<d4y2x2[3]<<std::endl;
    myfile.close();
}

void BoundaryDerivatives::outputZs(std::string filename){
    outputLeftRightZs(filename);
    outputTopBottomZs(filename); 
}

void BoundaryDerivatives::outputLeftRightZs(std::string filename){
    std::ofstream myfile;
    myfile.open("../outputs/tb_"+filename);
    myfile << "x,z_top,z_bottom"<<std::endl;
    for(int i = 0; i<m; i++){
        myfile<<x[i]<<","<<z_top[i]<<","<<z_bottom[i]<<std::endl;
    }
    myfile.close();
}

void BoundaryDerivatives::outputTopBottomZs(std::string filename){
    std::ofstream myfile;
    myfile.open("../outputs/lr_"+filename);
    myfile << "y,z_left,z_right"<<std::endl;
    for(int i = 0; i<n; i++){
        myfile<<y[i]<<","<<z_left[i]<<","<<z_right[i]<<std::endl;
    }
    myfile.close();
}