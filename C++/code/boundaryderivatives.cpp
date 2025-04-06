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
BoundaryDerivatives::BoundaryDerivatives(const std::vector<double>& Us, const std::vector<double>& Vs, const std::vector<std::vector<double>>& XFs) : us(Us), vs(Vs), xfs(XFs), m(us.size()), n(vs.size()), d2xfd2x_left(n), d2xfd2x_right(n), d2xfd2Q2_bottom(m), d2xfd2Q2_top(m), d4xfd2xd2Q2_vec(4), xfs_left(n), xfs_right(n), xfs_bottom(m), xfs_top(m) {
    // calculate boundary derivatives
    d2xfd2Q2();
    d2xfd2x();
    d4xfd2xd2Q2();
    // set them to struct
    SetBoundaryDerivatives();
}

// Second vs partial
void BoundaryDerivatives::d2xfd2Q2() {
    // Extrapolate to boundary 
    std::vector<double> tempVec1 = {vs[2],vs[1],vs[0]};
    std::vector<double> tempVec3 = {vs[n-3],vs[n-2],vs[n-1]};
    for (int i = 0; i < m; i++) {
        std::vector<double> tempVec2 = {xfs[2][i],xfs[1][i],xfs[0][i]};
        std::vector<double> tempVec4 = {xfs[n-3][i],xfs[n-2][i],xfs[n-1][i]};

        Extrapolate fit1(tempVec1, tempVec2);
        fit1.extrapolatePowerLaw(vs[0]);
        Extrapolate fit2(tempVec3, tempVec4);
        fit2.extrapolatePowerLaw(vs[n-1]);

        d2xfd2Q2_bottom[i] = fit1.solution;
        d2xfd2Q2_top[i] = fit2.solution;

        xfs_bottom[i] = fit1.powerLawxfs(vs[0]);
        xfs_top[i] = fit2.powerLawxfs(vs[n-1]);
    }
}

// Second us partial
void BoundaryDerivatives::d2xfd2x() {
    // Extrapolate to boundary 
    std::vector<double> tempVec1 = {us[1],us[0]};
    std::vector<double> tempVec3 = {us[m-2],us[m-1]};
    for (int j = 0; j < n; j++) {
        std::vector<double> tempVec2 = {xfs[j][1],xfs[j][0]}; 
        std::vector<double> tempVec4 = {xfs[j][m-3],xfs[j][m-2],xfs[j][m-1]};
        
        Extrapolate fit1(tempVec1, tempVec2);
        fit1.extrapolateTenPowerLaw(us[0]);
        Extrapolate fit2(tempVec3, tempVec4);
        fit2.extrapolateTenPowerLaw(us[m-1]);

        d2xfd2x_left[j] = fit1.solution;
        d2xfd2x_right[j] = fit2.solution;

        xfs_left[j] = fit1.tenPowerLawxfs(us[0]);
        xfs_right[j] = fit2.tenPowerLawxfs(us[m-1]);
    }
}

void BoundaryDerivatives::d4xfd2xd2Q2() {
    // Assume us and vs are independant, d4 is then product of d2s
    // corner derivatives go counter-clockwise: [0,0] [1,0] [1,1] [0,1]
    // Put in internal variable
    d4xfd2xd2Q2_vec[0] = d2xfd2x_left[0]*d2xfd2Q2_bottom[0];
    d4xfd2xd2Q2_vec[1] = d2xfd2x_right[0]*d2xfd2Q2_bottom[m-1];
    d4xfd2xd2Q2_vec[2] = d2xfd2x_right[n-1]*d2xfd2Q2_top[m-1];
    d4xfd2xd2Q2_vec[3] = d2xfd2x_left[n-1]*d2xfd2Q2_top[0];
}


// Sets the bounadry derivatives to the struct
void BoundaryDerivatives::SetBoundaryDerivatives(){
    boundary_derivatives.d2xfd2xs_left = d2xfd2x_left;
    boundary_derivatives.d2xfd2Q2s_top= d2xfd2Q2_top;
    boundary_derivatives.d2xfd2xs_right = d2xfd2x_right;
    boundary_derivatives.d2xfd2Q2s_bottom = d2xfd2Q2_bottom;
    boundary_derivatives.d4xfd2xd2Q2s_corners = d4xfd2xd2Q2_vec;
    
}

// Testing function that writes out the derivative values
void BoundaryDerivatives::outputDerivs(std::string filename){
    outputXDerivs(filename);
    outputQ2Derivs(filename); 
    outputCornerDerivs(filename);  
}

void BoundaryDerivatives::outputQ2Derivs(std::string filename){
    std::ofstream myfile;
    myfile.open("../outputs/y_"+filename);
    myfile << "us,d2xfQ2_top,d2xfQ2_bottom"<<std::endl;
    for(int i = 0; i<m; i++){
        myfile<<us[i]<<","<<d2xfd2Q2_top[i]<<","<<d2xfd2Q2_bottom[i]<<std::endl;
    }
    myfile.close();
}

void BoundaryDerivatives::outputXDerivs(std::string filename){
    std::ofstream myfile;
    myfile.open("../outputs/x_"+filename);
    myfile << "vs,d2xfd2x_left,d2xfd2x_right"<<std::endl;
    for(int i = 0; i<n; i++){
        myfile<<vs[i]<<","<<d2xfd2x_left[i]<<","<<d2xfd2x_right[i]<<std::endl;
    }
    myfile.close();
}

void BoundaryDerivatives::outputCornerDerivs(std::string filename){
    std::ofstream myfile;
    myfile.open("../outputs/corner_"+filename);
    myfile << "us,vs,d4xfd2xd2Q2"<<std::endl;
    myfile<<us[0]<<","<<vs[0]<<","<<d4xfd2xd2Q2_vec[0]<<std::endl;
    myfile<<us[m-1]<<","<<vs[0]<<","<<d4xfd2xd2Q2_vec[1]<<std::endl;
    myfile<<us[m-1]<<","<<vs[n-1]<<","<<d4xfd2xd2Q2_vec[2]<<std::endl;
    myfile<<us[0]<<","<<vs[n-1]<<","<<d4xfd2xd2Q2_vec[3]<<std::endl;
    myfile.close();
}

void BoundaryDerivatives::outputxfs(std::string filename){
    outputLeftRightxfs(filename);
    outputTopBottomxfs(filename); 
}

void BoundaryDerivatives::outputLeftRightxfs(std::string filename){
    std::ofstream myfile;
    myfile.open("../outputs/tb_"+filename);
    myfile << "us,xfs_top,xfs_bottom"<<std::endl;
    for(int i = 0; i<m; i++){
        myfile<<us[i]<<","<<xfs_top[i]<<","<<xfs_bottom[i]<<std::endl;
    }
    myfile.close();
}

void BoundaryDerivatives::outputTopBottomxfs(std::string filename){
    std::ofstream myfile;
    myfile.open("../outputs/lr_"+filename);
    myfile << "vs,xfs_left,xfs_right"<<std::endl;
    for(int i = 0; i<n; i++){
        myfile<<vs[i]<<","<<xfs_left[i]<<","<<xfs_right[i]<<std::endl;
    }
    myfile.close();
}