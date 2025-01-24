#include "../headers/bicubicspline.h"

// fast integer power function
inline constexpr double power(double base, int exponent) {
    if (exponent == 0) {
        return 1;
    } else if (exponent > 0) {
        return base * power(base, exponent - 1);
    } else {
        return 1 / power(base, -exponent);
    }
}

// fast square function
inline double square(double x){return x*x;}

// fast cube function
inline double cube(double x){return x*x*x;}

// constructor without edge derivatives -- set them all to 0
BicubicSpline::BicubicSpline(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<std::vector<double>>& Zs) : xs(Xs), ys(Ys), zs(Zs), m(Xs.size()), n(Ys.size()), d2x2s_left(n,0), d2x2s_right(n,0), d2y2s_bottom(m,0), d2y2s_top(m,0), d4x2y2s_corners(4,0), S(n-1, std::vector<Vector16>(m-1)){
    // check if list of function values has the appropriate length
    if ( zs.size() != n && zs[0].size() != m ){
        throw std::length_error( "Array of function values must have shape (m, n), where m and n are the grid dimensions along x and y, respectively!" );
    }
    // other checks / initialisations
    // calculate spline
    CalculateSpline();
}

// constructor with edge derivatives
BicubicSpline::BicubicSpline(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<std::vector<double>>& Zs, const Derivatives& ds) : xs(Xs), ys(Ys), zs(Zs), m(Xs.size()), n(Ys.size()), S(n-1, std::vector<Vector16>(m-1)){
    // check if list of function values has the appropriate length
    if ( zs.size() != n && zs[0].size() != m ){
        throw std::length_error( "Array of function values must have shape (m, n), where m and n are the grid dimensions along x and y, respectively!" );
    }
    // other checks / initialisations
    d2x2s_left = ds.d2x2s_left; 
    d2x2s_right = ds.d2x2s_right;
    d2y2s_bottom = ds.d2y2s_bottom;
    d2y2s_top = ds.d2y2s_top;
    d4x2y2s_corners = ds.d4x2y2s_corners;
    // calculate spline
    CalculateSpline();
}

void BicubicSpline::CalculateSpline(){
    // initialise matrices
    // matrix of x derivatives at knots
    std::vector<std::vector<double>> zs_x(n, std::vector<double>(m));

    // matrix of y derivatives at knots
    std::vector<std::vector<double>> zs_y(n, std::vector<double> (m));

    // matrix of xy derivatives at knots
    std::vector<std::vector<double>> zs_xy(n, std::vector<double> (m));

    // boundary vectors of xyy derivatives at top and bottom knots
    std::vector<double> zs_xyy_bottom(m);
    std::vector<double> zs_xyy_top(m);

    // construct 1D splines along x and y to find first 
    // and cross derivatives at spline knots
    // construct splines along boundary rows
    CubicSpline spline_bottom(xs, zs[0], d4x2y2s_corners[0], d4x2y2s_corners[1]);
    CubicSpline spline_top(xs, zs[n-1], d4x2y2s_corners[3], d4x2y2s_corners[2]);

    // get spline coefficients
    std::vector<std::vector<double>> S_bottom = spline_bottom.GetS();
    std::vector<std::vector<double>> S_top = spline_top.GetS();

    // calculate and store analytical x derivatives at boundary knots
    // the rest of the terms of the derivative are zero at the knots
    // due to (x - t) terms (see CubicSpline code)
    for (  int i = 0; i<zs_xyy_bottom.size()-1; i++){
        zs_xyy_bottom[i] = S_bottom[i][1];
    }
    zs_xyy_bottom[zs_xyy_bottom.size()-1] = S_bottom[S_bottom.size()-1][1] + 2*S_bottom[S_bottom.size()-1][2]*(xs[n-1] - xs[n-2]) + 3* S_bottom[S_bottom.size()-1][3]*power((xs[n-1] - xs[n-2]),2);
    for (  int i = 0; i<zs_xyy_top.size()-1; i++ ){
        zs_xyy_top[i] = S_top[i][1];
    }
    zs_xyy_top[zs_xyy_top.size()-1] = S_top[S_top.size()-1][1] + 2* S_top[S_top.size()-1][2]*(xs[n-1] - xs[n-2]) + 3*S_top[S_top.size()-1][3]*power((xs[n-1] - xs[n-2]),2);
    // loop over data rows
    for (  int i = 0; i<n; i++ ){
        
        // construct spline along ith row using the 
        // user input data (excluding the boundary rows)
        CubicSpline spline_x(xs, zs[i], d2x2s_left[i], d2x2s_right[i]);
        // get spline coefficients
        std::vector<std::vector<double>> S_x = spline_x.GetS();

        // calculate and store analytical x derivatives 
        // at knots along the ith row
        for (  int j = 0; j<zs_x[i].size()-1; j++ ){
            zs_x[i][j] = S_x[j][1];
        }
        zs_x[i][zs_x[i].size()-1] = S_x[S_x.size()-1][1] + 2*S_x[S_x.size()-1][2]*(xs[n-1] - xs[n-2]) + 3*S_x[S_x.size()-1][3]*power((xs[n-1] - xs[n-2]),2);
    }
    // loop over data columns
    for (  int i = 0; i < m; i++ ){
        // construct spline along jth column using the calculated 
        // spline derivatives (excluding boundary columns)
        std::vector<double> temp(zs_x.size()); 
        for (  int j = 0; j<zs_x.size(); j++ ){
            temp[j] = zs_x[j][i];
        }
        CubicSpline spline_xy(ys, temp, zs_xyy_bottom[i], zs_xyy_top[i]);
        // get spline coefficients
        std::vector<std::vector<double>> S_xy = spline_xy.GetS();

        // calculate and store analytical xy derivatives
        // at knots along the jth column
        for (  int j = 0; j<zs_xy.size()-1; j++ ){
            zs_xy[j][i] = S_xy[j][1];
        }
        zs_xy[zs_xy.size()-1][i] = S_xy[S_xy.size()-1][1] + 2*S_xy[S_xy.size()-1][2]*(ys[m-1] - ys[m-2]) + 3*S_xy[S_xy.size()-1][3]*power((ys[m-1] - ys[m-2]),2);
        
        // construct spline along jth column using the user input data
        // (excluding boundary columns since they have no further use)
        std::vector<double> temp2(zs.size());
        for (  int j = 0; j<zs.size(); j++ ){
            temp2[j] = zs[j][i];
        }
        CubicSpline spline_y(ys, temp2, d2y2s_bottom[i], d2y2s_top[i]);

        // get spline coefficients
        std::vector<std::vector<double>> S_y = spline_y.GetS();

        // calculate and store analytical xy derivatives
        // at knots along the jth column
        for (  int j = 0; j<zs_y.size()-1; j++ ){
            zs_y[j][i] = S_y[j][1];
        }
        zs_y[zs_y.size()-1][i] = S_y[S_y.size()-1][1] + 2*S_y[S_y.size()-1][2]*(ys[m-1] - ys[m-2]) + 3*S_y[S_y.size()-1][3]*power((ys[m-1] - ys[m-2]),2);
    }

    // fill knot values and derivatives into matrix
    // loop over the outer dimensions of the matrix
    for (  int j=0; j<n-1; j++ ){
        for (  int i=0; i<m-1; i++ ){
            // array of knot values and knot derivatives
            // at each iteration (grid cell)
            Vector16 Z;

            // function values
            Z[0] = zs[j][i];
            Z[1] = zs[j][i+1];
            Z[2] = zs[j+1][i];
            Z[3] = zs[j+1][i+1];

            // x derivatives
            Z[4] = zs_x[j][i];
            Z[5] = zs_x[j][i+1];
            Z[6] = zs_x[j+1][i];
            Z[7] = zs_x[j+1][i+1];

            // y derivatives
            Z[8] = zs_y[j][i];
            Z[9] = zs_y[j][i+1];
            Z[10] = zs_y[j+1][i];
            Z[11] = zs_y[j+1][i+1];

            // xy derivatives
            Z[12] = zs_xy[j][i];
            Z[13] = zs_xy[j][i+1];
            Z[14] = zs_xy[j+1][i];
            Z[15] = zs_xy[j+1][i+1];

            // matrix A such that A * S[j][i] = Z
            Matrix16 A;

            // loop over dimensions of A
            for (  int l=0; l<4; l++ ){
                for (  int k=0; k<4; k++ ){
                    // function value terms
                    A(0,4*l + k) = power(xs[i],k) * power(ys[j],l);
                    A(1,4*l + k) = power(xs[i+1],k) * power(ys[j],l);
                    A(2,4*l + k) = power(xs[i],k) * power(ys[j+1],l);
                    A(3,4*l + k) = power(xs[i+1],k) * power(ys[j+1],l);

                    //initialise remaining terms to zero
                    A(4,4*l + k) = A(5,4*l + k) = A(6,4*l + k) = A(7,4*l + k) = 0;
                    A(8,4*l + k) = A(9,4*l + k) = A(10,4*l + k) = A(11,4*l + k) = 0;
                    A(12,4*l + k) = A(13,4*l + k) = A(14,4*l + k) = A(15,4*l + k) = 0;

                    if ( k > 0 ){
                        // x derivative terms
                        A(4,4*l + k) = k * power(xs[i],k-1) * power(ys[j],l);
                        A(5,4*l + k) = k * power(xs[i+1],k-1) * power(ys[j],l);
                        A(6,4*l + k) = k * power(xs[i],k-1) * power(ys[j+1],l);
                        A(7,4*l + k) = k * power(xs[i+1],k-1) * power(ys[j+1],l);

                        if ( l > 0 ){
                            // y derivative terms
                            A(8,4*l + k) = l * power(xs[i],k) * power(ys[j],l-1);
                            A(9,4*l + k) = l * power(xs[i+1],k) * power(ys[j],l-1);
                            A(10,4*l + k) = l * power(xs[i],k) * power(ys[j+1],l-1);
                            A(11,4*l + k) = l * power(xs[i+1],k) * power(ys[j+1],l-1);

                            // xy derivative terms
                            A(12,4*l + k) = k * l * power(xs[i],k-1) * power(ys[j],l-1);
                            A(13,4*l + k) = k * l * power(xs[i+1],k-1) * power(ys[j],l-1);
                            A(14,4*l + k) = k * l * power(xs[i],k-1) * power(ys[j+1],l-1);
                            A(15,4*l + k) = k * l * power(xs[i+1],k-1) * power(ys[j+1],l-1);
                        }
                    } else if ( l > 0 ){
                        // remaining y derivative terms
                        A(8,4*l + k) = l * power(xs[i],k) * power(ys[j],l-1);
                        A(9,4*l + k) = l * power(xs[i+1],k) * power(ys[j],l-1);
                        A(10,4*l + k) = l * power(xs[i],k) * power(ys[j+1],l-1);
                        A(11,4*l + k) = l * power(xs[i+1],k) * power(ys[j+1],l-1);
                    }
                }
            }  

            // calculate the coefficient array for this iteration (grid cell)
            // meaning: S[j][i] = A^-1 * Z
            // S[j][i] = linalg::MatMul(A_inv, Z);
            S[j][i] = A.fullPivLu().solve(Z);
        }
    }
}

double BicubicSpline::evaluateSpline( double X, double Y ){
    /*
    Returns the value(s) of the interpolated function at x, y.

    params:  x:  x-coordinates(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
                y:  y-coordinate(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
    */
    double x = X;
    double y = Y;

    // initialise array storing function values
    double splines;

    // find indices before which each user input value should be
    // slotted into the arrays of knot coordinates such that the
    // arrays remain sorted; these indices (j, i) correspond to
    // the spline indices (i, j) at which to evaluate the spline
    // for each input point

    // relevant indices
    // np.searchsorted uses binary search
    // side='right': xs[i] <= x_eval < xs[i+1] as required
    // indices given are those of the knot values before which each
    // input value should be slotted in so need to be shifted back by one
    int ixs,iys;
    ixs = std::upper_bound(xs.begin(), xs.end(), x) - xs.begin() - 1;
    iys = std::upper_bound(ys.begin(), ys.end(), y) - ys.begin() - 1;
    // if final input value equal to (or greater than) the final knot 
    // value, rightmost value evaluates to final index of knot array + 1 
    // so it needs to be shifted back by one more
    if (x >= xs[xs.size()-1] ){ixs -= 1;}
    
    if (y >= ys[ys.size()-1] ){iys -= 1;}
    // evaluate spline at (j, i)-th coordinates and store in array
    splines = (           S[iys][ixs][0]                + S[iys][ixs][1]  * x
                        + S[iys][ixs][2]  * power(x,2)        + S[iys][ixs][3]  * power(x,3)
                        + S[iys][ixs][4]         * y    + S[iys][ixs][5]  * x    * y
                        + S[iys][ixs][6]  * power(x,2) * y    + S[iys][ixs][7]  * power(x,3) * y
                        + S[iys][ixs][8]         * power(y,2) + S[iys][ixs][9]  * x    * power(y,2)
                        + S[iys][ixs][10] * power(x,2) * power(y,2) + S[iys][ixs][11] * power(x,3) * power(y,2)
                        + S[iys][ixs][12]        * power(y,3) + S[iys][ixs][13] * x    * power(y,3)
                        + S[iys][ixs][14] * power(x,2) * power(y,3) + S[iys][ixs][15] * power(x,3) * power(y,3) );
    return splines;
}

std::vector<std::vector<double>> BicubicSpline::evaluateSpline( const std::vector<double>& x, const std::vector<double>& y ){
    /*
    Returns the value(s) of the interpolated function at x, y.

    params:  x:  x-coordinates(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
                y:  y-coordinate(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
    */
    int m,n;

    n = y.size();
    m = x.size();

    // initialise array storing function values
    std::vector<std::vector<double>> splines(n, std::vector<double>(m));

    // find indices before which each user input value should be
    // slotted into the arrays of knot coordinates such that the
    // arrays remain sorted; these indices (j, i) correspond to
    // the spline indices (i, j) at which to evaluate the spline
    // for each input point

    // relevant indices
    // np.searchsorted uses binary search
    // side='right': xs[i] <= x_eval < xs[i+1] as required
    // indices given are those of the knot values before which each
    // input value should be slotted in so need to be shifted back by one
    std::vector<double> ixs(m);
    std::vector<double> iys(n);
    for ( int i=0; i<m; i++ ){
        ixs[i] = std::upper_bound(xs.begin(), xs.end()-1, x[i]) - xs.begin() - 1;
    }
    for ( int i=0; i<n; i++ ){
        iys[i] = std::upper_bound(ys.begin(), ys.end()-1, y[i]) - ys.begin() - 1;
    }
        
    // if final input value equal to (or greater than) the final knot 
    // value, rightmost value evaluates to final index of knot array + 1 
    // so it needs to be shifted back by one more
    if (x[x.size()-1] >= xs[xs.size()-1] ){ixs[ixs.size()-1] -= 1;}
        
    if (y[y.size()-1] >= ys[ys.size()-1] ){iys[iys.size()-1] -= 1;}

    // loop over input coordinates
    for ( int i=0; i<n; i++ ){
        for ( int j=0; j<m; j++ ){
            // evaluate spline at (j, i)-th coordinates and store in array
            splines[i][j] = (   S[iys[i]][ixs[j]][0]                + S[iys[i]][ixs[j]][1]  * x[j]
                                + S[iys[i]][ixs[j]][2]  * power(x[j],2)        + S[iys[i]][ixs[j]][3]  * power(x[j],3)
                                + S[iys[i]][ixs[j]][4]         * y[i]    + S[iys[i]][ixs[j]][5]  * x[j]    * y[i]
                                + S[iys[i]][ixs[j]][6]  * power(x[j],2) * y[i]    + S[iys[i]][ixs[j]][7]  * power(x[j],3) * y[i]
                                + S[iys[i]][ixs[j]][8]         * power(y[i],2) + S[iys[i]][ixs[j]][9]  * x[j]    * power(y[i],2)
                                + S[iys[i]][ixs[j]][10] * power(x[j],2) * power(y[i],2) + S[iys[i]][ixs[j]][11] * power(x[j],3) * power(y[i],2)
                                + S[iys[i]][ixs[j]][12]        * power(y[i],3) + S[iys[i]][ixs[j]][13] * x[j]    * power(y[i],3)
                                + S[iys[i]][ixs[j]][14] * power(x[j],2) * power(y[i],3) + S[iys[i]][ixs[j]][15] * power(x[j],3) * power(y[i],3) );
        }
    }
        
    return splines;
}

void BicubicSpline::calculateDerivs(const std::vector<double>& X, const std::vector<double>& Y){
    initDerivs(X, Y);

    std::vector<int> idxs(dm);
    std::vector<int> idys(dn);

    for ( int i=0; i<dm; i++ ){
        idxs[i] = std::upper_bound(xs.begin(), xs.end()-1, dxs[i]) - xs.begin() - 1;
    }

    for ( int i=0; i<dn; i++ ){
        idys[i] = std::upper_bound(ys.begin(), ys.end()-1, dys[i]) - ys.begin() - 1;
    }

    if (X[X.size()-1] >= xs[xs.size()-1] ){idxs[idxs.size()-1] -= 1;}
        
    if (Y[Y.size()-1] >= ys[ys.size()-1] ){idys[idys.size()-1] -= 1;}

    // Calculate each derivatvive
    d1X(dxs,dys,idxs,idys);
    d2X(dxs,dys,idxs,idys);
    d3X(dxs,dys,idxs,idys);
    d1Y(dxs,dys,idxs,idys);
    d2Y(dxs,dys,idxs,idys);
    d3Y(dxs,dys,idxs,idys);
}

void BicubicSpline::initDerivs(const std::vector<double>& X, const std::vector<double>& Y){
    // set internal variables
    // grid
    dxs = X;
    dys = Y;

    // grid dimensions
    dm = dxs.size();
    dn = dys.size();

    // initialise the derivative vectors
    d1x.reserve(dn);
    d2x.reserve(dn);
    d3x.reserve(dn);
    d1y.reserve(dn);
    d2y.reserve(dn);
    d3y.reserve(dn);
    for (int i = 0; i<dn; i++){
        d1x.push_back(std::vector<double>(dm));
        d2x.push_back(std::vector<double>(dm));
        d3x.push_back(std::vector<double>(dm));
        d1y.push_back(std::vector<double>(dm));
        d2y.push_back(std::vector<double>(dm));
        d3y.push_back(std::vector<double>(dm));
        }
}

void BicubicSpline::outputDerivs(const std::string& filename){
    std::ofstream myfile;
    myfile.open("../outputs/"+filename);
    myfile << "x,y,d1x,d1y,d2x,d2y,d3x,d3y"<<std::endl;

    for ( int i = 0; i<dn; i++ ){
        for( int j = 0; j<dm; j++ ){
            myfile<<dxs[j]<<","<<dys[i]<<","<<d1x[i][j]<<","<<d1y[i][j]<<","<<d2x[i][j]<<","<<d2y[i][j]<<","<<d3x[i][j]<<","<<d3y[i][j]<<std::endl;
        }
    }
    myfile.close();
}

// Derivative calculations
void BicubicSpline::d1X(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<int>& iXs, const std::vector<int>& iYs){
    int iy,ix;
    double y, x, result;

    for (int j = 0; j<dn; j++){
        iy = iYs[j];
        y = Ys[j];
        for (int i = 0; i<dm; i++){
            ix = iXs[i];
            x = Xs[i];

            result = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result += k*S[iy][ix][k + 4*l]*power(x,k-1)*power(y,l);
                }
            }

            d1x[j][i] = result / (std::pow(10,x)*std::log(10));
        }
    }
}

void BicubicSpline::d2X(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<int>& iXs, const std::vector<int>& iYs){
    int iy,ix;
    double y, x, result1, result2;

    for (int j = 0; j<dn; j++){
        iy = iYs[j];
        y = Ys[j];
        for (int i = 0; i<dm; i++){
            ix = iXs[i];
            x = Xs[i];

            result1 = 0;
            result2 = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result1 += k*(k-1)*S[iy][ix][k + 4*l]*power(x,k-2)*power(y,l);
                    result2 += k*S[iy][ix][k + 4*l]*power(x,k-1)*power(y,l);
                }
            }
            
            d2x[j][i] = (result1 / (std::pow(10,2*x)*power(std::log(10),2))) - (result2 / (std::pow(10,2*x)*std::log(10)));
        }
    }
}

void BicubicSpline::d3X(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<int>& iXs, const std::vector<int>& iYs){
    int iy,ix;
    double y, x, result1, result2, result3;

    for (int j = 0; j<dn; j++){
        iy = iYs[j];
        y = Ys[j];
        for (int i = 0; i<dm; i++){
            ix = iXs[i];
            x = Xs[i];
            
            result1 = 0;
            result2 = 0;
            result3 = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result1 += k*(k-1)*(k-2)*S[iy][ix][k + 4*l]*power(x,k-3)*power(y,l);
                    result2 += k*(k-1)*S[iy][ix][k + 4*l]*power(x,k-2)*power(y,l);
                    result3 += k*S[iy][ix][k + 4*l]*power(x,k-1)*power(y,l);
                }
            }
            
            d3x[j][i] = (result1 / (std::pow(10,3*x)*power(std::log(10),3))) - (3*result2 / (std::pow(10,3*x)*power(std::log(10),2))) + (2*result3 / (std::pow(10,3*x)*std::log(10)));
        }
    }
}

void BicubicSpline::d1Y(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<int>& iXs, const std::vector<int>& iYs){
    int iy,ix;
    double y, x, result;

    for (int j = 0; j<dn; j++){
        iy = iYs[j];
        y = Ys[j];
        for (int i = 0; i<dm; i++){
            ix = iXs[i];
            x = Xs[i];

            result = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result += l*S[iy][ix][l*4 + k]*power(x,k)*power(y,l-1);
                }
            }
            
            d1y[j][i] = result / (std::pow(10,y)*std::log(10));
        }
    }
}

void BicubicSpline::d2Y(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<int>& iXs, const std::vector<int>& iYs){
    int iy,ix;
    double y, x, result1, result2;

    for (int j = 0; j<dn; j++){
        iy = iYs[j];
        y = Ys[j];
        for (int i = 0; i<dm; i++){
            ix = iXs[i];
            x = Xs[i];
            
            result1 = 0;
            result2 = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result1 += l*(l-1)*S[iy][ix][k + 4*l]*power(x,k)*power(y,l-2);
                    result2 += l*S[iy][ix][k + 4*l]*power(x,k)*power(y,l-1);
                }
            }
            
            d2y[j][i] = (result1 / (std::pow(10,2*y)*power(std::log(10),2))) - (result2 / (std::pow(10,2*y)*std::log(10)));
        }
    }
}

void BicubicSpline::d3Y(const std::vector<double>& Xs, const std::vector<double>& Ys, const std::vector<int>& iXs, const std::vector<int>& iYs){
    int iy,ix;
    double y, x, result1, result2, result3;

    for (int j = 0; j<dn; j++){
        iy = iYs[j];
        y = Ys[j];
        for (int i = 0; i<dm; i++){
            ix = iXs[i];
            x = Xs[i];

            result1 = 0;
            result2 = 0;
            result3 = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result1 += l*(l-1)*(l-2)*S[iy][ix][k + 4*l]*power(x,k)*power(y,l-3);
                    result2 += l*(l-1)*S[iy][ix][k + 4*l]*power(x,k)*power(y,l-2);
                    result3 += l*S[iy][ix][k + 4*l]*power(x,k)*power(y,l-1);
                }
            }
            
            d3y[j][i] = (result1 / (std::pow(10,3*y)*power(std::log(10),3))) - (3*result2 / (std::pow(10,3*y)*power(std::log(10),2))) + (2*result3 / (std::pow(10,3*y)*std::log(10)));
        }
    }
}