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
BicubicSpline::BicubicSpline(const std::vector<double>& _us, const std::vector<double>& _vs, const std::vector<std::vector<double>>& _xfs) : us(_us), vs(_vs), xfs(_xfs), m(_us.size()), n(_vs.size()), d2xfd2xs_left(n,0), d2xfd2xs_right(n,0), d2xfd2Q2s_bottom(m,0), d2xfd2Q2s_top(m,0), d4xfd2xd2Q2s_corners(4,0), S((n-1)*(m-1)*16){
    // check if list of function values has the appropriate length
    if ( xfs.size() != n || xfs[0].size() != m ){
        throw std::length_error( "Array of function values must have shape (n, m), where m and n are the grid dimensions along x and y, respectively!" );
    }
    // other checks / initialisations
    // calculate spline
    CalculateSpline();
}

// constructor with edge derivatives
BicubicSpline::BicubicSpline(const std::vector<double>& _us, const std::vector<double>& _vs, const std::vector<std::vector<double>>& _xfs, const Derivatives& ds) : us(_us), vs(_vs), xfs(_xfs), m(_us.size()), n(_vs.size()), S((n-1)*(m-1)*16){
    // check if list of function values has the appropriate length
    if ( xfs.size() != n || xfs[0].size() != m ){
        throw std::length_error( "Array of function values must have shape (m, n), where m and n are the grid dimensions along x and y, respectively!" );
    }
    // other checks / initialisations
    d2xfd2xs_left = ds.d2xfd2xs_left; 
    d2xfd2xs_right = ds.d2xfd2xs_right;
    d2xfd2Q2s_bottom = ds.d2xfd2Q2s_bottom;
    d2xfd2Q2s_top = ds.d2xfd2Q2s_top;
    d4xfd2xd2Q2s_corners = ds.d4xfd2xd2Q2s_corners;
    // calculate spline
    CalculateSpline();
}

void BicubicSpline::CalculateSpline(){
    // initialise matrices
    // matrix of x derivatives at knots
    std::vector<std::vector<double>> xfs_x(n, std::vector<double>(m));

    // matrix of y derivatives at knots
    std::vector<std::vector<double>> xfs_y(n, std::vector<double> (m));

    // matrix of xy derivatives at knots
    std::vector<std::vector<double>> xfs_xy(n, std::vector<double> (m));

    // boundary vectors of xyy derivatives at top and bottom knots
    std::vector<double> xfs_xyy_bottom(m);
    std::vector<double> xfs_xyy_top(m);

    // construct 1D splines along x and y to find first 
    // and cross derivatives at spline knots
    // construct splines along boundary rows
    CubicSpline spline_bottom(us, xfs[0], d4xfd2xd2Q2s_corners[0], d4xfd2xd2Q2s_corners[1]);
    CubicSpline spline_top(us, xfs[n-1], d4xfd2xd2Q2s_corners[3], d4xfd2xd2Q2s_corners[2]);

    // get spline coefficients
    std::vector<std::vector<double>> S_bottom = spline_bottom.GetS();
    std::vector<std::vector<double>> S_top = spline_top.GetS();

    // calculate and store analytical x derivatives at boundary knots
    // the rest of the terms of the derivative are zero at the knots
    // due to (x - t) terms (see CubicSpline code)
    for (  int i = 0; i<xfs_xyy_bottom.size()-1; i++){
        xfs_xyy_bottom[i] = S_bottom[i][1];
    }
    xfs_xyy_bottom[xfs_xyy_bottom.size()-1] = S_bottom[S_bottom.size()-1][1] + 2*S_bottom[S_bottom.size()-1][2]*(us[n-1] - us[n-2]) + 3* S_bottom[S_bottom.size()-1][3]*power((us[n-1] - us[n-2]),2);
    for (  int i = 0; i<xfs_xyy_top.size()-1; i++ ){
        xfs_xyy_top[i] = S_top[i][1];
    }
    xfs_xyy_top[xfs_xyy_top.size()-1] = S_top[S_top.size()-1][1] + 2* S_top[S_top.size()-1][2]*(us[n-1] - us[n-2]) + 3*S_top[S_top.size()-1][3]*power((us[n-1] - us[n-2]),2);
    // loop over data rows
    for (  int i = 0; i<n; i++ ){
        
        // construct spline along ith row using the 
        // user input data (excluding the boundary rows)
        CubicSpline spline_x(us, xfs[i], d2xfd2xs_left[i], d2xfd2xs_right[i]);
        // get spline coefficients
        std::vector<std::vector<double>> S_x = spline_x.GetS();

        // calculate and store analytical x derivatives 
        // at knots along the ith row
        for (  int j = 0; j<xfs_x[i].size()-1; j++ ){
            xfs_x[i][j] = S_x[j][1];
        }
        xfs_x[i][xfs_x[i].size()-1] = S_x[S_x.size()-1][1] + 2*S_x[S_x.size()-1][2]*(us[n-1] - us[n-2]) + 3*S_x[S_x.size()-1][3]*power((us[n-1] - us[n-2]),2);
    }
    // loop over data columns
    for (  int i = 0; i < m; i++ ){
        // construct spline along jth column using the calculated 
        // spline derivatives (excluding boundary columns)
        std::vector<double> temp(xfs_x.size()); 
        for (  int j = 0; j<xfs_x.size(); j++ ){
            temp[j] = xfs_x[j][i];
        }
        CubicSpline spline_xy(vs, temp, xfs_xyy_bottom[i], xfs_xyy_top[i]);
        // get spline coefficients
        std::vector<std::vector<double>> S_xy = spline_xy.GetS();

        // calculate and store analytical xy derivatives
        // at knots along the jth column
        for (  int j = 0; j<xfs_xy.size()-1; j++ ){
            xfs_xy[j][i] = S_xy[j][1];
        }
        xfs_xy[xfs_xy.size()-1][i] = S_xy[S_xy.size()-1][1] + 2*S_xy[S_xy.size()-1][2]*(vs[m-1] - vs[m-2]) + 3*S_xy[S_xy.size()-1][3]*power((vs[m-1] - vs[m-2]),2);
        
        // construct spline along jth column using the user input data
        // (excluding boundary columns since they have no further use)
        std::vector<double> temp2(xfs.size());
        for (  int j = 0; j<xfs.size(); j++ ){
            temp2[j] = xfs[j][i];
        }
        CubicSpline spline_y(vs, temp2, d2xfd2Q2s_bottom[i], d2xfd2Q2s_top[i]);

        // get spline coefficients
        std::vector<std::vector<double>> S_y = spline_y.GetS();

        // calculate and store analytical xy derivatives
        // at knots along the jth column
        for (  int j = 0; j<xfs_y.size()-1; j++ ){
            xfs_y[j][i] = S_y[j][1];
        }
        xfs_y[xfs_y.size()-1][i] = S_y[S_y.size()-1][1] + 2*S_y[S_y.size()-1][2]*(vs[m-1] - vs[m-2]) + 3*S_y[S_y.size()-1][3]*power((vs[m-1] - vs[m-2]),2);
    }

    // fill knot values and derivatives into matrix
    // loop over the outer dimensions of the matrix
    for (  int j=0; j<n-1; j++ ){
        for (  int i=0; i<m-1; i++ ){
            // array of knot values and knot derivatives
            // at each iteration (grid cell)
            Vector16 Z;

            // function values
            Z[0] = xfs[j][i];
            Z[1] = xfs[j][i+1];
            Z[2] = xfs[j+1][i];
            Z[3] = xfs[j+1][i+1];

            // x derivatives
            Z[4] = xfs_x[j][i];
            Z[5] = xfs_x[j][i+1];
            Z[6] = xfs_x[j+1][i];
            Z[7] = xfs_x[j+1][i+1];

            // y derivatives
            Z[8] = xfs_y[j][i];
            Z[9] = xfs_y[j][i+1];
            Z[10] = xfs_y[j+1][i];
            Z[11] = xfs_y[j+1][i+1];

            // xy derivatives
            Z[12] = xfs_xy[j][i];
            Z[13] = xfs_xy[j][i+1];
            Z[14] = xfs_xy[j+1][i];
            Z[15] = xfs_xy[j+1][i+1];

            // matrix A such that A * S[j][i] = Z
            Matrix16 A;

            // loop over dimensions of A
            for (  int l=0; l<4; l++ ){
                for (  int k=0; k<4; k++ ){
                    // function value terms
                    A(0,4*l + k) = power(us[i],k) * power(vs[j],l);
                    A(1,4*l + k) = power(us[i+1],k) * power(vs[j],l);
                    A(2,4*l + k) = power(us[i],k) * power(vs[j+1],l);
                    A(3,4*l + k) = power(us[i+1],k) * power(vs[j+1],l);

                    //initialise remaining terms to zero
                    A(4,4*l + k) = A(5,4*l + k) = A(6,4*l + k) = A(7,4*l + k) = 0;
                    A(8,4*l + k) = A(9,4*l + k) = A(10,4*l + k) = A(11,4*l + k) = 0;
                    A(12,4*l + k) = A(13,4*l + k) = A(14,4*l + k) = A(15,4*l + k) = 0;

                    if ( k > 0 ){
                        // x derivative terms
                        A(4,4*l + k) = k * power(us[i],k-1) * power(vs[j],l);
                        A(5,4*l + k) = k * power(us[i+1],k-1) * power(vs[j],l);
                        A(6,4*l + k) = k * power(us[i],k-1) * power(vs[j+1],l);
                        A(7,4*l + k) = k * power(us[i+1],k-1) * power(vs[j+1],l);

                        if ( l > 0 ){
                            // y derivative terms
                            A(8,4*l + k) = l * power(us[i],k) * power(vs[j],l-1);
                            A(9,4*l + k) = l * power(us[i+1],k) * power(vs[j],l-1);
                            A(10,4*l + k) = l * power(us[i],k) * power(vs[j+1],l-1);
                            A(11,4*l + k) = l * power(us[i+1],k) * power(vs[j+1],l-1);

                            // xy derivative terms
                            A(12,4*l + k) = k * l * power(us[i],k-1) * power(vs[j],l-1);
                            A(13,4*l + k) = k * l * power(us[i+1],k-1) * power(vs[j],l-1);
                            A(14,4*l + k) = k * l * power(us[i],k-1) * power(vs[j+1],l-1);
                            A(15,4*l + k) = k * l * power(us[i+1],k-1) * power(vs[j+1],l-1);
                        }
                    } else if ( l > 0 ){
                        // remaining y derivative terms
                        A(8,4*l + k) = l * power(us[i],k) * power(vs[j],l-1);
                        A(9,4*l + k) = l * power(us[i+1],k) * power(vs[j],l-1);
                        A(10,4*l + k) = l * power(us[i],k) * power(vs[j+1],l-1);
                        A(11,4*l + k) = l * power(us[i+1],k) * power(vs[j+1],l-1);
                    }
                }
            }  

            // calculate the coefficient array for this iteration (grid cell)
            // meaning: S[j][i] = A^-1 * Z
            // S[j][i] = linalg::MatMul(A_inv, Z);
            Vector16 temp = A.completeOrthogonalDecomposition().pseudoInverse()*Z;
            for (size_t k = 0; k<16; ++k){
                S[j*(m-1)*16 + i*16 + k] = temp[k];
            }    
        }
    }
}

double BicubicSpline::evaluateSpline( double X, double Y ) const {
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
    double splines, result;

    // find indices before which each user input value should be
    // slotted into the arravs of knot coordinates such that the
    // arravs remain sorted; these indices (j, i) correspond to
    // the spline indices (i, j) at which to evaluate the spline
    // for each input point

    // relevant indices
    // np.searchsorted uses binary search
    // side='right': us[i] <= x_eval < us[i+1] as required
    // indices given are those of the knot values before which each
    // input value should be slotted in so need to be shifted back by one
    int ix,iy;
    ix = std::upper_bound(us.begin(), us.end(), x) - us.begin() - 1;
    iy = std::upper_bound(vs.begin(), vs.end(), y) - vs.begin() - 1;
    // if final input value equal to (or greater than) the final knot 
    // value, rightmost value evaluates to final index of knot array + 1 
    // so it needs to be shifted back by one more
    if (x >= us[us.size()-1] ){ix -= 1;}
    
    if (y >= vs[vs.size()-1] ){iy -= 1;}

    result = 0;
    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += S[iy*(m-1)*16 + ix*16 + 4*l + k]*power(x,k)*power(y,l);
        }
    }
    
    // evaluate spline at (j, i)-th coordinates and store in array
    splines = result;
    return splines;
}

double BicubicSpline::evaluateSpline( double x, double y, size_t ix, size_t iy ) const {
    // initialise array storing function values
    const double x2 = x*x;
    const double x3 = x2*x;
    const double y2 = y*y;
    const double y3 = y2*y;
    const uint irow_col = (iy*(m - 1)+ ix)*16;


    const double result = S[irow_col] + S[irow_col + 4] * y + S[irow_col + 8] * y2 + S[irow_col + 12] * y3 
           + S[irow_col + 1] * x + S[irow_col + 5] * x * y + S[irow_col + 9] * x * y2 + S[irow_col + 13] * x * y3 
           + S[irow_col + 2] * x2 + S[irow_col + 6] * x2 * y + S[irow_col + 10] * x2 * y2 + S[irow_col + 14] * x2 * y3 
           + S[irow_col + 3] * x3 + S[irow_col + 7] * x3 * y + S[irow_col + 11] * x3 * y2 + S[irow_col + 15] * x3 * y3;
    
    return result;
}

std::vector<std::vector<double>> BicubicSpline::evaluateSpline( const std::vector<double>& x, const std::vector<double>& y ) const {
    /*
    Returns the value(s) of the interpolated function at x, y.

    params:  x:  x-coordinates(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
                y:  y-coordinate(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
    */
    int m, n, iy, ix;
    double y_val, x_val, result;

    n = y.size();
    m = x.size();

    // initialise array storing function values
    std::vector<std::vector<double>> splines(n, std::vector<double>(m));

    // find indices before which each user input value should be
    // slotted into the arravs of knot coordinates such that the
    // arravs remain sorted; these indices (j, i) correspond to
    // the spline indices (i, j) at which to evaluate the spline
    // for each input point

    // relevant indices
    // np.searchsorted uses binary search
    // side='right': us[i] <= x_eval < us[i+1] as required
    // indices given are those of the knot values before which each
    // input value should be slotted in so need to be shifted back by one
    std::vector<double> ius(m);
    std::vector<double> ivs(n);
    for ( int i=0; i<m; i++ ){
        ius[i] = std::upper_bound(us.begin(), us.end()-1, x[i]) - us.begin() - 1;
    }
    for ( int i=0; i<n; i++ ){
        ivs[i] = std::upper_bound(vs.begin(), vs.end()-1, y[i]) - vs.begin() - 1;
    }
        
    // if final input value equal to (or greater than) the final knot 
    // value, rightmost value evaluates to final index of knot array + 1 
    // so it needs to be shifted back by one more
    if (x[x.size()-1] >= us[us.size()-1] ){ius[ius.size()-1] -= 1;}
        
    if (y[y.size()-1] >= vs[vs.size()-1] ){ivs[ivs.size()-1] -= 1;}

    // loop over input coordinates
    for ( int i=0; i<n; i++ ){
        iy = ivs[i];
        y_val = y[i];
        for ( int j=0; j<m; j++ ){
            ix = ius[j];
            x_val = x[j];

            result = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result += S[iy*(BicubicSpline::m - 1)*16 + ix*16 + 4*l + k]*power(x_val,k)*power(y_val,l);
                }
            }
            // evaluate spline at (j, i)-th coordinates and store in array
            splines[i][j] = result;
        }
    }
        
    return splines;
}

void BicubicSpline::calculateDerivs(const std::vector<double>& X, const std::vector<double>& Y){
    initDerivs(X, Y);

    std::vector<int> idus(dm);
    std::vector<int> idvs(dn);

    for ( int i=0; i<dm; i++ ){
        idus[i] = std::upper_bound(us.begin(), us.end()-1, dus[i]) - us.begin() - 1;
    }

    for ( int i=0; i<dn; i++ ){
        idvs[i] = std::upper_bound(vs.begin(), vs.end()-1, dvs[i]) - vs.begin() - 1;
    }

    if (X[X.size()-1] >= us[us.size()-1] ){idus[idus.size()-1] -= 1;}
        
    if (Y[Y.size()-1] >= vs[vs.size()-1] ){idvs[idvs.size()-1] -= 1;}

    // Calculate each derivative
    dxfdx(dus,dvs,idus,idvs);
    d2xfd2x(dus,dvs,idus,idvs);
    d3xfd3x(dus,dvs,idus,idvs);
    dxfdQ2(dus,dvs,idus,idvs);
    d2xfd2Q2(dus,dvs,idus,idvs);
    d3xfd3Q2(dus,dvs,idus,idvs);
}

std::tuple<double,double,double,double> BicubicSpline::evaluateDerivs(double du, double dv, size_t idu, size_t idv) const{
    double _dxfdx,_dxfdQ2,_d2xfd2x,_d2xfd2Q2;
    // Calculate each derivative
    _dxfdx = dxfdx(du,dv,idu,idv);
    _d2xfd2x = d2xfd2x(du,dv,idu,idv);
    _dxfdQ2 = dxfdQ2(du,dv,idu,idv);
    _d2xfd2Q2 = d2xfd2Q2(du,dv,idu,idv);
    return  std::make_tuple(_dxfdx, _dxfdQ2, _d2xfd2x, _d2xfd2Q2);
}

std::tuple<double,double,double,double> BicubicSpline::evaluateLogDerivs(double du, double dv, size_t idu, size_t idv) const{
    double _dxfdu,_dxfdv,_d2xfd2u,_d2xfd2v;
    // Calculate each derivative
    _dxfdu = dxfdu(du,dv,idu,idv);
    _d2xfd2u = d2xfd2u(du,dv,idu,idv);
    _dxfdv = dxfdv(du,dv,idu,idv);
    _d2xfd2v = d2xfd2v(du,dv,idu,idv);
    return  std::make_tuple(_dxfdu, _dxfdv, _d2xfd2u, _d2xfd2v);
}

void BicubicSpline::initDerivs(const std::vector<double>& X, const std::vector<double>& Y){
    // set internal variables
    // grid
    dus = X;
    dvs = Y;

    // grid dimensions
    dm = dus.size();
    dn = dvs.size();

    // initialise the derivative vectors
    dxfdx_vec.reserve(dn);
    d2xfd2x_vec.reserve(dn);
    d3xfd3x_vec.reserve(dn);
    dxfdQ2_vec.reserve(dn);
    d2xfd2Q2_vec.reserve(dn);
    d3xfd3Q2_vec.reserve(dn);
    for (int i = 0; i<dn; i++){
        dxfdx_vec.push_back(std::vector<double>(dm));
        d2xfd2x_vec.push_back(std::vector<double>(dm));
        d3xfd3x_vec.push_back(std::vector<double>(dm));
        dxfdQ2_vec.push_back(std::vector<double>(dm));
        d2xfd2Q2_vec.push_back(std::vector<double>(dm));
        d3xfd3Q2_vec.push_back(std::vector<double>(dm));
        }
}

void BicubicSpline::outputDerivs(const std::string& filename) const {
    std::ofstream myfile;
    myfile.open("../outputs/"+filename);
    myfile << "x,y,dxfdx,dxfdQ2,d2xfd2x,d2xfd2Q2,d3xfd3x,d3xfd3Q2_vec"<<std::endl;

    for ( int i = 0; i<dn; i++ ){
        for( int j = 0; j<dm; j++ ){
            myfile<<dus[j]<<","<<dvs[i]<<","<<dxfdx_vec[i][j]<<","<<dxfdQ2_vec[i][j]<<","<<d2xfd2x_vec[i][j]<<","<<d2xfd2Q2_vec[i][j]<<","<<d3xfd3x_vec[i][j]<<","<<d3xfd3Q2_vec[i][j]<<std::endl;
        }
    }
    myfile.close();
}

// Derivative calculations
void BicubicSpline::dxfdx(const std::vector<double>& _us, const std::vector<double>& _vs, const std::vector<int>& i_us, const std::vector<int>& i_vs){
    int iy,ix;
    double y, x, result;

    for (int j = 0; j<dn; j++){
        iy = i_vs[j];
        y = _vs[j];
        for (int i = 0; i<dm; i++){
            ix = i_us[i];
            x = _us[i];

            result = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result += k*S[iy*(BicubicSpline::m-1)*16 + ix*16 + 4*l + k]*power(x,k-1)*power(y,l);
                }
            }

            dxfdx_vec[j][i] = result / (std::pow(10,x)*std::log(10));
        }
    }
}

void BicubicSpline::d2xfd2x(const std::vector<double>& _us, const std::vector<double>& _vs, const std::vector<int>& i_us, const std::vector<int>& i_vs){
    int iy,ix;
    double y, x, result1, result2;

    for (int j = 0; j<dn; j++){
        iy = i_vs[j];
        y = _vs[j];
        for (int i = 0; i<dm; i++){
            ix = i_us[i];
            x = _us[i];

            result1 = 0;
            result2 = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result1 += k*(k-1)*S[iy*(BicubicSpline::m-1)*16 + ix*16 + 4*l + k]*power(x,k-2)*power(y,l);
                    result2 += k*S[iy*(BicubicSpline::m-1)*16 + ix*16 + 4*l + k]*power(x,k-1)*power(y,l);
                }
            }
            
            d2xfd2x_vec[j][i] = (result1 / (std::pow(10,2*x)*power(std::log(10),2))) - (result2 / (std::pow(10,2*x)*std::log(10)));
        }
    }
}

void BicubicSpline::d3xfd3x(const std::vector<double>& _us, const std::vector<double>& _vs, const std::vector<int>& i_us, const std::vector<int>& i_vs){
    int iy,ix;
    double y, x, result1, result2, result3;

    for (int j = 0; j<dn; j++){
        iy = i_vs[j];
        y = _vs[j];
        for (int i = 0; i<dm; i++){
            ix = i_us[i];
            x = _us[i];
            
            result1 = 0;
            result2 = 0;
            result3 = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result1 += k*(k-1)*(k-2)*S[iy*(BicubicSpline::m-1)*16 + ix*16 + 4*l + k]*power(x,k-3)*power(y,l);
                    result2 += k*(k-1)*S[iy*(BicubicSpline::m-1)*16 + ix*16 + 4*l + k]*power(x,k-2)*power(y,l);
                    result3 += k*S[iy*(BicubicSpline::m-1)*16 + ix*16 + 4*l + k]*power(x,k-1)*power(y,l);
                }
            }
            
            d3xfd3x_vec[j][i] = (result1 / (std::pow(10,3*x)*power(std::log(10),3))) - (3*result2 / (std::pow(10,3*x)*power(std::log(10),2))) + (2*result3 / (std::pow(10,3*x)*std::log(10)));
        }
    }
}

void BicubicSpline::dxfdQ2(const std::vector<double>& _us, const std::vector<double>& _vs, const std::vector<int>& i_us, const std::vector<int>& i_vs){
    int iy,ix;
    double y, x, result;

    for (int j = 0; j<dn; j++){
        iy = i_vs[j];
        y = _vs[j];
        for (int i = 0; i<dm; i++){
            ix = i_us[i];
            x = _us[i];

            result = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result += l*S[iy*(BicubicSpline::m-1)*16 + ix*16 + 4*l + k]*power(x,k)*power(y,l-1);
                }
            }
            
            dxfdQ2_vec[j][i] = result / (std::pow(10,y)*std::log(10));
        }
    }
}

void BicubicSpline::d2xfd2Q2(const std::vector<double>& _us, const std::vector<double>& _vs, const std::vector<int>& i_us, const std::vector<int>& i_vs){
    int iy,ix;
    double y, x, result1, result2;

    for (int j = 0; j<dn; j++){
        iy = i_vs[j];
        y = _vs[j];
        for (int i = 0; i<dm; i++){
            ix = i_us[i];
            x = _us[i];
            
            result1 = 0;
            result2 = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result1 += l*(l-1)*S[iy*(BicubicSpline::m-1)*16 + ix*16 + 4*l + k]*power(x,k)*power(y,l-2);
                    result2 += l*S[iy*(BicubicSpline::m-1)*16 + ix*16 + 4*l + k]*power(x,k)*power(y,l-1);
                }
            }
            
            d2xfd2Q2_vec[j][i] = (result1 / (std::pow(10,2*y)*power(std::log(10),2))) - (result2 / (std::pow(10,2*y)*std::log(10)));
        }
    }
}

void BicubicSpline::d3xfd3Q2(const std::vector<double>& _us, const std::vector<double>& _vs, const std::vector<int>& i_us, const std::vector<int>& i_vs){
    int iy,ix;
    double y, x, result1, result2, result3;

    for (int j = 0; j<dn; j++){
        iy = i_vs[j];
        y = _vs[j];
        for (int i = 0; i<dm; i++){
            ix = i_us[i];
            x = _us[i];

            result1 = 0;
            result2 = 0;
            result3 = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result1 += l*(l-1)*(l-2)*S[iy*(BicubicSpline::m-1)*16 + ix*16 + 4*l + k]*power(x,k)*power(y,l-3);
                    result2 += l*(l-1)*S[iy*(BicubicSpline::m-1)*16 + ix*16 + 4*l + k]*power(x,k)*power(y,l-2);
                    result3 += l*S[iy*(BicubicSpline::m-1)*16 + ix*16 + 4*l + k]*power(x,k)*power(y,l-1);
                }
            }
            
            d3xfd3Q2_vec[j][i] = (result1 / (std::pow(10,3*y)*power(std::log(10),3))) - (3*result2 / (std::pow(10,3*y)*power(std::log(10),2))) + (2*result3 / (std::pow(10,3*y)*std::log(10)));
        }
    }
}

double BicubicSpline::dxfdx(double x, double y, size_t ix, size_t iy) const{
    double result = 0;
    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += k*S[iy*(m-1)*16 + ix*16 + 4*l + k]*power(x,k-1)*power(y,l);
        }
    }

    return(result / std::exp(x));
}

double BicubicSpline::d2xfd2x(double x, double y, size_t ix, size_t iy) const{
    double result1 = 0;
    double result2 = 0;

    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result1 += k*(k-1)*S[iy*(m-1)*16 + ix*16 + 4*l + k]*power(x,k-2)*power(y,l);
            result2 += k*S[iy*(m-1)*16 + ix*16 + 4*l + k]*power(x,k-1)*power(y,l);
        }
    }
    
    return((result1 - result2) / std::exp(2*x));
}

double BicubicSpline::dxfdQ2(double x, double y, size_t ix, size_t iy) const{
    double result = 0;

    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += l*S[iy*(m-1)*16 + ix*16 + 4*l + k]*power(x,k)*power(y,l-1);
        }
    }
    
    return(result / std::exp(y));
}

double BicubicSpline::d2xfd2Q2(double x, double y, size_t ix, size_t iy) const{
    double result1 = 0;
    double result2 = 0;

    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result1 += l*(l-1)*S[iy*(m-1)*16 + ix*16 + 4*l + k]*power(x,k)*power(y,l-2);
            result2 += l*S[iy*(m-1)*16 + ix*16 + 4*l + k]*power(x,k)*power(y,l-1);
        }
    }
    
    return((result1 - result2)/ std::exp(2*y));
}

double BicubicSpline::dxfdu(double x, double y, size_t ix, size_t iy) const{
    double result = 0;
    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += k*S[iy*(m-1)*16 + ix*16 + 4*l + k]*power(x,k-1)*power(y,l);
        }
    }

    return result;
}

double BicubicSpline::d2xfd2u(double x, double y, size_t ix, size_t iy) const{
    double result = 0;

    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += k*(k-1)*S[iy*(m-1)*16 + ix*16 + 4*l + k]*power(x,k-2)*power(y,l);
        }
    }
    
    return result;
}

double BicubicSpline::dxfdv(double x, double y, size_t ix, size_t iy) const{
    double result = 0;

    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += l*S[iy*(m-1)*16 + ix*16 + 4*l + k]*power(x,k)*power(y,l-1);
        }
    }
    
    return result;
}

double BicubicSpline::d2xfd2v(double x, double y, size_t ix, size_t iy) const{
    double result = 0;

    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += l*(l-1)*S[iy*(m-1)*16 + ix*16 + 4*l + k]*power(x,k)*power(y,l-2);
        }
    }
    
    return result;
}