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
BicubicSpline::BicubicSpline(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<std::vector<double>>& xfs) : us(us), vs(vs), xfs(xfs), m(us.size()), n(vs.size()), d2xfd2xs_left(n,0), d2xfd2xs_right(n,0), d2xfd2Q2s_bottom(m,0), d2xfd2Q2s_top(m,0), d4xfd2xd2Q2s_corners(4,0), S((n-1)*(m-1)*16){
    // check if list of function values has the appropriate length
    if ( xfs.size() != n || xfs[0].size() != m ){
        throw std::length_error( "Array of function values must have shape (n, m), where m and n are the grid dimensions along u and v, respectively!" );
    }
    // other checks / initialisations
    // calculate spline
    CalculateSpline();
}

// constructor with edge derivatives
BicubicSpline::BicubicSpline(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<std::vector<double>>& xfs, const Derivatives& ds) : us(us), vs(vs), xfs(xfs), m(us.size()), n(vs.size()), S((n-1)*(m-1)*16){
    // check if list of function values has the appropriate length
    if ( xfs.size() != n || xfs[0].size() != m ){
        throw std::length_error( "Array of function values must have shape (m, n), where m and n are the grid dimensions along u and v, respectively!" );
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
    // matrix of u derivatives at knots
    std::vector<std::vector<double>> dxfsdx(n, std::vector<double>(m));

    // matrix of v derivatives at knots
    std::vector<std::vector<double>> dxfsdQ2(n, std::vector<double> (m));

    // matrix of xy derivatives at knots
    std::vector<std::vector<double>> d2xfsdxdQ2(n, std::vector<double> (m));

    // boundary vectors of xyy derivatives at top and bottom knots
    std::vector<double> xfs_xyy_bottom(m);
    std::vector<double> xfs_xyy_top(m);

    // construct 1D splines along u and v to find first 
    // and cross derivatives at spline knots
    // construct splines along boundary rows
    CubicSpline spline_bottom(us, xfs[0], d4xfd2xd2Q2s_corners[0], d4xfd2xd2Q2s_corners[1]);
    CubicSpline spline_top(us, xfs[n-1], d4xfd2xd2Q2s_corners[3], d4xfd2xd2Q2s_corners[2]);

    // get spline coefficients
    std::vector<std::vector<double>> S_bottom = spline_bottom.GetS();
    std::vector<std::vector<double>> S_top = spline_top.GetS();

    // calculate and store analytical u derivatives at boundary knots
    // the rest of the terms of the derivative are zero at the knots
    // due to (u - t) terms (see CubicSpline code)
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
        CubicSpline spline_u(us, xfs[i], d2xfd2xs_left[i], d2xfd2xs_right[i]);
        // get spline coefficients
        std::vector<std::vector<double>> S_u = spline_u.GetS();

        // calculate and store analytical u derivatives 
        // at knots along the ith row
        for (  int j = 0; j<dxfsdx[i].size()-1; j++ ){
            dxfsdx[i][j] = S_u[j][1];
        }
        dxfsdx[i][dxfsdx[i].size()-1] = S_u[S_u.size()-1][1] + 2*S_u[S_u.size()-1][2]*(us[n-1] - us[n-2]) + 3*S_u[S_u.size()-1][3]*power((us[n-1] - us[n-2]),2);
    }
    // loop over data columns
    for (  int i = 0; i < m; i++ ){
        // construct spline along jth column using the calculated 
        // spline derivatives (excluding boundary columns)
        std::vector<double> temp(dxfsdx.size()); 
        for (  int j = 0; j<dxfsdx.size(); j++ ){
            temp[j] = dxfsdx[j][i];
        }
        CubicSpline spline_uv(vs, temp, xfs_xyy_bottom[i], xfs_xyy_top[i]);
        // get spline coefficients
        std::vector<std::vector<double>> S_uv = spline_uv.GetS();

        // calculate and store analytical xy derivatives
        // at knots along the jth column
        for (  int j = 0; j<d2xfsdxdQ2.size()-1; j++ ){
            d2xfsdxdQ2[j][i] = S_uv[j][1];
        }
        d2xfsdxdQ2[d2xfsdxdQ2.size()-1][i] = S_uv[S_uv.size()-1][1] + 2*S_uv[S_uv.size()-1][2]*(vs[m-1] - vs[m-2]) + 3*S_uv[S_uv.size()-1][3]*power((vs[m-1] - vs[m-2]),2);
        
        // construct spline along jth column using the user input data
        // (excluding boundary columns since they have no further use)
        std::vector<double> temp2(xfs.size());
        for (  int j = 0; j<xfs.size(); j++ ){
            temp2[j] = xfs[j][i];
        }
        CubicSpline spline_v(vs, temp2, d2xfd2Q2s_bottom[i], d2xfd2Q2s_top[i]);

        // get spline coefficients
        std::vector<std::vector<double>> S_v = spline_v.GetS();

        // calculate and store analytical xy derivatives
        // at knots along the jth column
        for (  int j = 0; j<dxfsdQ2.size()-1; j++ ){
            dxfsdQ2[j][i] = S_v[j][1];
        }
        dxfsdQ2[dxfsdQ2.size()-1][i] = S_v[S_v.size()-1][1] + 2*S_v[S_v.size()-1][2]*(vs[m-1] - vs[m-2]) + 3*S_v[S_v.size()-1][3]*power((vs[m-1] - vs[m-2]),2);
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

            // u derivatives
            Z[4] = dxfsdx[j][i];
            Z[5] = dxfsdx[j][i+1];
            Z[6] = dxfsdx[j+1][i];
            Z[7] = dxfsdx[j+1][i+1];

            // v derivatives
            Z[8] = dxfsdQ2[j][i];
            Z[9] = dxfsdQ2[j][i+1];
            Z[10] = dxfsdQ2[j+1][i];
            Z[11] = dxfsdQ2[j+1][i+1];

            // xy derivatives
            Z[12] = d2xfsdxdQ2[j][i];
            Z[13] = d2xfsdxdQ2[j][i+1];
            Z[14] = d2xfsdxdQ2[j+1][i];
            Z[15] = d2xfsdxdQ2[j+1][i+1];

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
                        // u derivative terms
                        A(4,4*l + k) = k * power(us[i],k-1) * power(vs[j],l);
                        A(5,4*l + k) = k * power(us[i+1],k-1) * power(vs[j],l);
                        A(6,4*l + k) = k * power(us[i],k-1) * power(vs[j+1],l);
                        A(7,4*l + k) = k * power(us[i+1],k-1) * power(vs[j+1],l);

                        if ( l > 0 ){
                            // v derivative terms
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
                        // remaining v derivative terms
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

double BicubicSpline::evaluateSpline( double u, double v ) const {
    /*
    Returns the value(s) of the interpolated function at u, v.

    params:  u:  u-coordinates(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
                v:  v-coordinate(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
    */
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
    int iu,iv;
    iu = std::upper_bound(us.begin(), us.end(), u) - us.begin() - 1;
    iv = std::upper_bound(vs.begin(), vs.end(), v) - vs.begin() - 1;
    // if final input value equal to (or greater than) the final knot 
    // value, rightmost value evaluates to final index of knot array + 1 
    // so it needs to be shifted back by one more
    if (u >= us[us.size()-1] ){iu -= 1;}
    
    if (v >= vs[vs.size()-1] ){iv -= 1;}

    result = 0;
    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += S[iv*(m-1)*16 + iu*16 + 4*l + k]*power(u,k)*power(v,l);
        }
    }
    
    // evaluate spline at (j, i)-th coordinates and store in array
    splines = result;
    return splines;
}

double BicubicSpline::evaluateSpline( double u, double v, size_t iu, size_t iv ) const {
    // initialise array storing function values
    const double u2 = u*u;
    const double u3 = u2*u;
    const double v2 = v*v;
    const double v3 = v2*v;
    const uint irow_col = (iv*(m - 1)+ iu)*16;


    const double result = S[irow_col] + S[irow_col + 4] * v + S[irow_col + 8] * v2 + S[irow_col + 12] * v3 
           + S[irow_col + 1] * u + S[irow_col + 5] * u * v + S[irow_col + 9] * u * v2 + S[irow_col + 13] * u * v3 
           + S[irow_col + 2] * u2 + S[irow_col + 6] * u2 * v + S[irow_col + 10] * u2 * v2 + S[irow_col + 14] * u2 * v3 
           + S[irow_col + 3] * u3 + S[irow_col + 7] * u3 * v + S[irow_col + 11] * u3 * v2 + S[irow_col + 15] * u3 * v3;
    
    return result;
}

std::vector<std::vector<double>> BicubicSpline::evaluateSpline( const std::vector<double>& u, const std::vector<double>& v ) const {
    /*
    Returns the value(s) of the interpolated function at u, v.

    params:  u:  u-coordinates(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
                v:  v-coordinate(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
    */
    int m, n, iv, iu;
    double y_val, x_val, result;

    n = v.size();
    m = u.size();

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
        ius[i] = std::upper_bound(us.begin(), us.end()-1, u[i]) - us.begin() - 1;
    }
    for ( int i=0; i<n; i++ ){
        ivs[i] = std::upper_bound(vs.begin(), vs.end()-1, v[i]) - vs.begin() - 1;
    }
        
    // if final input value equal to (or greater than) the final knot 
    // value, rightmost value evaluates to final index of knot array + 1 
    // so it needs to be shifted back by one more
    if (u[u.size()-1] >= us[us.size()-1] ){ius[ius.size()-1] -= 1;}
        
    if (v[v.size()-1] >= vs[vs.size()-1] ){ivs[ivs.size()-1] -= 1;}

    // loop over input coordinates
    for ( int i=0; i<n; i++ ){
        iv = ivs[i];
        y_val = v[i];
        for ( int j=0; j<m; j++ ){
            iu = ius[j];
            x_val = u[j];

            result = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result += S[iv*(BicubicSpline::m - 1)*16 + iu*16 + 4*l + k]*power(x_val,k)*power(y_val,l);
                }
            }
            // evaluate spline at (j, i)-th coordinates and store in array
            splines[i][j] = result;
        }
    }
        
    return splines;
}

void BicubicSpline::calculateDerivs(const std::vector<double>& u, const std::vector<double>& v){
    initDerivs(u, v);

    std::vector<int> idus(dm);
    std::vector<int> idvs(dn);

    for ( int i=0; i<dm; i++ ){
        idus[i] = std::upper_bound(us.begin(), us.end()-1, dus[i]) - us.begin() - 1;
    }

    for ( int i=0; i<dn; i++ ){
        idvs[i] = std::upper_bound(vs.begin(), vs.end()-1, dvs[i]) - vs.begin() - 1;
    }

    if (u[u.size()-1] >= us[us.size()-1] ){idus[idus.size()-1] -= 1;}
        
    if (v[v.size()-1] >= vs[vs.size()-1] ){idvs[idvs.size()-1] -= 1;}

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

void BicubicSpline::initDerivs(const std::vector<double>& du, const std::vector<double>& dv){
    // set internal variables
    // grid
    dus = du;
    dvs = dv;

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
    myfile << "u,v,dxfdx,dxfdQ2,d2xfd2x,d2xfd2Q2,d3xfd3x,d3xfd3Q2_vec"<<std::endl;

    for ( int i = 0; i<dn; i++ ){
        for( int j = 0; j<dm; j++ ){
            myfile<<dus[j]<<","<<dvs[i]<<","<<dxfdx_vec[i][j]<<","<<dxfdQ2_vec[i][j]<<","<<d2xfd2x_vec[i][j]<<","<<d2xfd2Q2_vec[i][j]<<","<<d3xfd3x_vec[i][j]<<","<<d3xfd3Q2_vec[i][j]<<std::endl;
        }
    }
    myfile.close();
}

// Derivative calculations
void BicubicSpline::dxfdx(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<int>& ius, const std::vector<int>& ivs){
    int iv,iu;
    double v, u, result;

    for (int j = 0; j<dn; j++){
        iv = ivs[j];
        v = vs[j];
        for (int i = 0; i<dm; i++){
            iu = ius[i];
            u = us[i];

            result = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result += k*S[iv*(BicubicSpline::m-1)*16 + iu*16 + 4*l + k]*power(u,k-1)*power(v,l);
                }
            }

            dxfdx_vec[j][i] = result / (std::pow(10,u)*std::log(10));
        }
    }
}

void BicubicSpline::d2xfd2x(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<int>& ius, const std::vector<int>& ivs){
    int iv,iu;
    double v, u, result1, result2;

    for (int j = 0; j<dn; j++){
        iv = ivs[j];
        v = vs[j];
        for (int i = 0; i<dm; i++){
            iu = ius[i];
            u = us[i];

            result1 = 0;
            result2 = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result1 += k*(k-1)*S[iv*(BicubicSpline::m-1)*16 + iu*16 + 4*l + k]*power(u,k-2)*power(v,l);
                    result2 += k*S[iv*(BicubicSpline::m-1)*16 + iu*16 + 4*l + k]*power(u,k-1)*power(v,l);
                }
            }
            
            d2xfd2x_vec[j][i] = (result1 / (std::pow(10,2*u)*power(std::log(10),2))) - (result2 / (std::pow(10,2*u)*std::log(10)));
        }
    }
}

void BicubicSpline::d3xfd3x(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<int>& ius, const std::vector<int>& ivs){
    int iv,iu;
    double v, u, result1, result2, result3;

    for (int j = 0; j<dn; j++){
        iv = ivs[j];
        v = vs[j];
        for (int i = 0; i<dm; i++){
            iu = ius[i];
            u = us[i];
            
            result1 = 0;
            result2 = 0;
            result3 = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result1 += k*(k-1)*(k-2)*S[iv*(BicubicSpline::m-1)*16 + iu*16 + 4*l + k]*power(u,k-3)*power(v,l);
                    result2 += k*(k-1)*S[iv*(BicubicSpline::m-1)*16 + iu*16 + 4*l + k]*power(u,k-2)*power(v,l);
                    result3 += k*S[iv*(BicubicSpline::m-1)*16 + iu*16 + 4*l + k]*power(u,k-1)*power(v,l);
                }
            }
            
            d3xfd3x_vec[j][i] = (result1 / (std::pow(10,3*u)*power(std::log(10),3))) - (3*result2 / (std::pow(10,3*u)*power(std::log(10),2))) + (2*result3 / (std::pow(10,3*u)*std::log(10)));
        }
    }
}

void BicubicSpline::dxfdQ2(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<int>& ius, const std::vector<int>& ivs){
    int iv,iu;
    double v, u, result;

    for (int j = 0; j<dn; j++){
        iv = ivs[j];
        v = vs[j];
        for (int i = 0; i<dm; i++){
            iu = ius[i];
            u = us[i];

            result = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result += l*S[iv*(BicubicSpline::m-1)*16 + iu*16 + 4*l + k]*power(u,k)*power(v,l-1);
                }
            }
            
            dxfdQ2_vec[j][i] = result / (std::pow(10,v)*std::log(10));
        }
    }
}

void BicubicSpline::d2xfd2Q2(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<int>& ius, const std::vector<int>& ivs){
    int iv,iu;
    double v, u, result1, result2;

    for (int j = 0; j<dn; j++){
        iv = ivs[j];
        v = vs[j];
        for (int i = 0; i<dm; i++){
            iu = ius[i];
            u = us[i];
            
            result1 = 0;
            result2 = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result1 += l*(l-1)*S[iv*(BicubicSpline::m-1)*16 + iu*16 + 4*l + k]*power(u,k)*power(v,l-2);
                    result2 += l*S[iv*(BicubicSpline::m-1)*16 + iu*16 + 4*l + k]*power(u,k)*power(v,l-1);
                }
            }
            
            d2xfd2Q2_vec[j][i] = (result1 / (std::pow(10,2*v)*power(std::log(10),2))) - (result2 / (std::pow(10,2*v)*std::log(10)));
        }
    }
}

void BicubicSpline::d3xfd3Q2(const std::vector<double>& us, const std::vector<double>& vs, const std::vector<int>& ius, const std::vector<int>& ivs){
    int iv,iu;
    double v, u, result1, result2, result3;

    for (int j = 0; j<dn; j++){
        iv = ivs[j];
        v = vs[j];
        for (int i = 0; i<dm; i++){
            iu = ius[i];
            u = us[i];

            result1 = 0;
            result2 = 0;
            result3 = 0;
            for (int k=0; k<4; k++){
                for (int l=0; l<4; l++){
                    result1 += l*(l-1)*(l-2)*S[iv*(BicubicSpline::m-1)*16 + iu*16 + 4*l + k]*power(u,k)*power(v,l-3);
                    result2 += l*(l-1)*S[iv*(BicubicSpline::m-1)*16 + iu*16 + 4*l + k]*power(u,k)*power(v,l-2);
                    result3 += l*S[iv*(BicubicSpline::m-1)*16 + iu*16 + 4*l + k]*power(u,k)*power(v,l-1);
                }
            }
            
            d3xfd3Q2_vec[j][i] = (result1 / (std::pow(10,3*v)*power(std::log(10),3))) - (3*result2 / (std::pow(10,3*v)*power(std::log(10),2))) + (2*result3 / (std::pow(10,3*v)*std::log(10)));
        }
    }
}

double BicubicSpline::dxfdx(double u, double v, size_t iu, size_t iv) const{
    double result = 0;
    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += k*S[iv*(m-1)*16 + iu*16 + 4*l + k]*power(u,k-1)*power(v,l);
        }
    }

    return(result / std::exp(u));
}

double BicubicSpline::d2xfd2x(double u, double v, size_t iu, size_t iv) const{
    double result1 = 0;
    double result2 = 0;

    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result1 += k*(k-1)*S[iv*(m-1)*16 + iu*16 + 4*l + k]*power(u,k-2)*power(v,l);
            result2 += k*S[iv*(m-1)*16 + iu*16 + 4*l + k]*power(u,k-1)*power(v,l);
        }
    }
    
    return((result1 - result2) / std::exp(2*u));
}

double BicubicSpline::dxfdQ2(double u, double v, size_t iu, size_t iv) const{
    double result = 0;

    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += l*S[iv*(m-1)*16 + iu*16 + 4*l + k]*power(u,k)*power(v,l-1);
        }
    }
    
    return(result / std::exp(v));
}

double BicubicSpline::d2xfd2Q2(double u, double v, size_t iu, size_t iv) const{
    double result1 = 0;
    double result2 = 0;

    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result1 += l*(l-1)*S[iv*(m-1)*16 + iu*16 + 4*l + k]*power(u,k)*power(v,l-2);
            result2 += l*S[iv*(m-1)*16 + iu*16 + 4*l + k]*power(u,k)*power(v,l-1);
        }
    }
    
    return((result1 - result2)/ std::exp(2*v));
}

double BicubicSpline::dxfdu(double u, double v, size_t iu, size_t iv) const{
    double result = 0;
    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += k*S[iv*(m-1)*16 + iu*16 + 4*l + k]*power(u,k-1)*power(v,l);
        }
    }

    return result;
}

double BicubicSpline::d2xfd2u(double u, double v, size_t iu, size_t iv) const{
    double result = 0;

    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += k*(k-1)*S[iv*(m-1)*16 + iu*16 + 4*l + k]*power(u,k-2)*power(v,l);
        }
    }
    
    return result;
}

double BicubicSpline::dxfdv(double u, double v, size_t iu, size_t iv) const{
    double result = 0;

    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += l*S[iv*(m-1)*16 + iu*16 + 4*l + k]*power(u,k)*power(v,l-1);
        }
    }
    
    return result;
}

double BicubicSpline::d2xfd2v(double u, double v, size_t iu, size_t iv) const{
    double result = 0;

    for (int k=0; k<4; k++){
        for (int l=0; l<4; l++){
            result += l*(l-1)*S[iv*(m-1)*16 + iu*16 + 4*l + k]*power(u,k)*power(v,l-2);
        }
    }
    
    return result;
}