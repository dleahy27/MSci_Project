#include "../headers/finitederivatives.h"

inline double square(double x){return x*x;}

inline constexpr double power(double base, int exponent) {
    if (exponent == 0) {
        return 1;
    } else if (exponent > 0) {
        return base * power(base, exponent - 1);
    } else {
        return std::pow(base,exponent);
    }
}

FiniteDerivatives::FiniteDerivatives(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<std::vector<double>>& Z) : x(X), y(Y), z(Z), m(x.size()), n(y.size()), d1x(n, std::vector<double>(m)), d2x(n, std::vector<double>(m)), d1y(n, std::vector<double>(m)), d2y(n, std::vector<double>(m)), d2_left(n), d2_right(n), d2_bottom(m), d2_top(m), d4y2x2(4) {
    finiteD1();
    finiteD2();
    finiteBoundary();
    SetBoundaryDerivatives();
}

void FiniteDerivatives::finiteD1(){
    // Loop through each point in the matrix
    partialDerivative1X();
    partialDerivative1Y();
}

void FiniteDerivatives::finiteD2(){
    // Loop through each point in the matrix
    partialDerivative2X();
    partialDerivative2Y();
    partialDerivative2X2Y();
}

// Function to calculate the partial derivative with respect to x at each z point
void FiniteDerivatives::partialDerivative1X() {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (j == 0) {
                // Forward difference at the left edge
                double hx_forward = x[j + 1] - x[j];
                d1x[i][j] = (z[i][j + 1] - z[i][j]) / hx_forward;
            } else if (j == m - 1) {
                // Backward difference at the right edge
                double hx_backward = x[j] - x[j - 1];
                d1x[i][j] = (z[i][j] - z[i][j - 1]) / hx_backward;
            } else {
                // Central difference for internal points
                double hx_forward = x[j + 1] - x[j];
                double hx_backward = x[j] - x[j - 1];
                d1x[i][j] = (square(hx_backward)*z[i][j + 1] + (square(hx_forward) - square(hx_backward))*z[i][j] - square(hx_forward)*z[i][j - 1]) / (hx_forward*hx_backward*(hx_forward + hx_backward));
            }
        }
    }
}

// Function to calculate the partial derivative with respect to y at each z point
void FiniteDerivatives::partialDerivative1Y() {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (i == 0) {
                // Forward difference at the top edge
                double hy_forward = y[i + 1] - y[i];
                d1y[i][j] = (z[i + 1][j] - z[i][j]) / hy_forward;
            } else if (i == n - 1) {
                // Backward difference at the bottom edge
                double hy_backward = y[i] - y[i - 1];
                d1y[i][j] = (z[i][j] - z[i - 1][j]) / hy_backward;
            } else {
                // Central difference for internal points
                double hy_forward = y[i + 1] - y[i];
                double hy_backward = y[i] - y[i - 1];
                d1y[i][j] = (square(hy_backward)*z[i + 1][j] + (square(hy_forward) - square(hy_backward))*z[i][j] - square(hy_forward)*z[i - 1][j]) / (hy_forward*hy_backward*(hy_forward + hy_backward));
            }
        }
    }
}

void FiniteDerivatives::partialDerivative2Y() {
    for (int i = 1; i < n-1; ++i) {
        for (int j = 0; j < m; ++j) {
            if (i == 1) {
                // Forward difference at the top edge
                double hy_forward = y[i + 1] - y[i];
                d2y[i][j] = (d1y[i + 1][j] - d1y[i][j]) / hy_forward;
            } else if (i == n - 2) {
                // Backward difference at the bottom edge
                double hy_backward = y[i] - y[i - 1];
                d2y[i][j] = (d1y[i][j] - d1y[i - 1][j]) / hy_backward;
            } else {
                // Central difference for internal points
                double hy_forward = y[i + 1] - y[i];
                double hy_backward = y[i] - y[i - 1];
                d2y[i][j] = (square(hy_backward)*d1y[i + 1][j] + (square(hy_forward) - square(hy_backward))*d1y[i][j] - square(hy_forward)*d1y[i - 1][j]) / (hy_forward*hy_backward*(hy_forward + hy_backward));
            }
        }    
    }

    for (int i = 1; i < m-1; i++) {
        std::vector<double> tempVec1 = {y[3],y[2],y[1]}, tempVec2 = {d2y[3][i],d2y[2][i],d2y[1][i]};
        std::vector<double> tempVec3 = {y[n-4],y[n-3],y[n-2]}, tempVec4 = {d2y[n-4][i],d2y[n-3][i],d2y[n-2][i]};

        Extrapolate fit1(y[0], tempVec1, tempVec2);
        Extrapolate fit2(y[n-1], tempVec3, tempVec4);

        d2y[0][i] = fit1.solution;
        d2y[n-1][i] = fit2.solution;
    }
}

void FiniteDerivatives::partialDerivative2X() {
    for (int i = 0; i < n; i++) {
        for (int j = 1; j < m-1; j++) {
            if (j == 1) {
                // Forward difference at the left edge
                double hx_forward = x[j + 1] - x[j];
                d2x[i][j] = (d1x[i][j + 1] - d1x[i][j]) / hx_forward;
            } else if (j == n - 2) {
                // Backward difference at the right edge
                double hx_backward = x[j] - x[j - 1];
                d2x[i][j] = (d1x[i][j] - d1x[i][j - 1]) / hx_backward;
            } else {
                // Central difference for internal points
                double hx_forward = x[j + 1] - x[j];
                double hx_backward = x[j] - x[j - 1];
                d2x[i][j] = (square(hx_backward)*d1x[i][j + 1] + (square(hx_forward) - square(hx_backward))*d1x[i][j] - square(hx_forward)*d1x[i][j - 1]) / (hx_forward*hx_backward*(hx_forward + hx_backward));
            }
        }
    }

    std::vector<double> tempVec1 = {x[3],x[2],x[1]};
    std::vector<double> tempVec3 = {x[m-4],x[m-3],x[m-2]};
    for (int j = 1; j < n-1; j++) {
        std::vector<double> tempVec2 = {d2x[j][3],d2x[j][2],d2x[j][1]};
        std::vector<double> tempVec4 = {d2x[j][m-4],d2x[j][m-3],d2x[j][m-2]};
        
        Extrapolate fit1(x[0], tempVec1, tempVec2);
        Extrapolate fit2(x[m-1], tempVec3, tempVec4);

        d2x[j][0] = fit1.solution;
        d2x[j][m-1] = fit2.solution;
    }
}

void FiniteDerivatives::partialDerivative2X2Y() {
    // Assume x and y are independant, d4 is then product of d2s
    // corner derivatives go counter-clockwise: [0,0] [1,0] [1,1] [0,1]
    std::vector<double> XS_begin = {x[3],x[2],x[1]};
    std::vector<double> YS_begin = {y[3],y[2],y[1]};
    std::vector<double> XS_end = {x[m-4],x[m-3],x[m-2]};
    std::vector<double> YS_end = {y[n-4],y[n-3],y[n-2]};

    // Extrapolate to the corners
    // xs
    Extrapolate xfit1(y[0], YS_begin, {d2x[3][0], d2x[2][0], d2x[1][0]});
    d2x[0][0] = xfit1.solution;

    Extrapolate xfit2(y[n-1], YS_end, {d2x[n-4][0], d2x[n-3][0], d2x[n-2][0]});
    d2x[n-1][0] = xfit2.solution;

    Extrapolate xfit3(y[n-1], YS_end, {d2x[n-4][m-1], d2x[n-3][m-1], d2x[n-2][m-1]});
    d2x[n-1][m-1] = xfit3.solution;

    Extrapolate xfit4(y[0], YS_begin, {d2x[3][m-1], d2x[2][m-1], d2x[1][m-1]});
    d2x[0][m-1] = xfit4.solution;

    // ys
    Extrapolate yfit1(x[0], XS_begin, {d2y[0][3], d2y[0][2], d2y[0][1]});
    d2y[0][0] = yfit1.solution;

    Extrapolate yfit2(x[0], XS_begin, {d2y[n-1][3], d2y[n-1][2], d2y[n-1][1]});
    d2y[n-1][0] = yfit2.solution;

    Extrapolate yfit3(x[m-1], XS_end, {d2y[n-1][m-4], d2y[n-1][m-3], d2y[n-1][m-2]});
    d2y[n-1][m-1] = yfit3.solution;

    Extrapolate yfit4(x[m-1], XS_end, {d2y[0][m-4], d2y[0][m-3], d2y[0][m-2]});
    d2y[0][m-1] = yfit4.solution;
    
    // Put in internal variable
    d4y2x2[0] = d2x[0][0]*d2y[0][0];
    d4y2x2[1] = d2x[0][m-1]*d2y[0][m-1];
    d4y2x2[2] = d2x[n-1][m-1]*d2y[n-1][m-1];
    d4y2x2[3] = d2x[n-1][0]*d2y[n-1][0];
}

void FiniteDerivatives::d2Left(){
    for ( int i = 0; i < n ; i++){
        d2_left[i] = d2x[i][0];
    }
}

void FiniteDerivatives::d2Right(){
    for ( int i = 0; i < n; i++){
        d2_right[i] = d2x[i][m-1];
    }
}

void FiniteDerivatives::d2Top(){
    for ( int i = 0; i < m; i++){
        d2_top[i] = d2y[n-1][i];
    }
}

void FiniteDerivatives::d2Bottom(){
    for ( int i = 0; i < m; i++){
        d2_bottom[i] = d2y[0][i];
    }
}

void FiniteDerivatives::finiteBoundary(){
    d2Left();
    d2Right();
    d2Top();
    d2Bottom();
}

void FiniteDerivatives::SetBoundaryDerivatives(){
    boundary_derivatives.d2x2s_left = d2_left;
    boundary_derivatives.d2y2s_top = d2_top;
    boundary_derivatives.d2x2s_right = d2_right;
    boundary_derivatives.d2y2s_bottom = d2_bottom;
    boundary_derivatives.d4x2y2s_corners = d4y2x2;
    
}

void FiniteDerivatives::outputDerivs(std::string filename){
    
    std::ofstream myfile;
    myfile.open("../outputs/"+filename);
    myfile << "x,y,d1x,d1y,d2x,d2y"<<std::endl;

    for ( int i = 0; i<n; i++ ){
        for( int j = 0; j<m; j++ ){
            myfile<<x[j]<<","<<y[i]<<","<<d1x[i][j]<<","<<d1y[i][j]<<","<<d2x[i][j]<<","<<d2y[i][j]<<std::endl;
        }
    }
    myfile.close();
}
