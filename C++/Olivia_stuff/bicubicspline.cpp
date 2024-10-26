#include "bicubicspline.h"

BicubicSpline::BicubicSpline(std::vector<double>Xs, std::vector<double>Ys, std::vector<std::vector<double>>Zs) : xs(Xs), ys(Ys), zs(Zs), m(Xs.size()), n(Ys.size()), d2x2s_left(n,0), d2x2s_right(n,0), d2y2s_bottom(m,0), d2y2s_top(m,0), d4x2y2s_corners(4,0), _S(n-1, std::vector<std::vector<double>>(m-1, std::vector<double>(16, 0.0))){
    bool verbose = true;
    if ( verbose ){
        std::cout << "In constructor:" << std::endl;
    }

    if ( verbose ){
        std::cout << "Performing checks..." << std::endl; 
    }

    // check if list of function values has the appropriate length
    if ( zs.size() != n && zs[0].size() != m ){
        throw std::length_error( "Array of function values must have shape (m, n), where m and n are the grid dimensions along x and y, respectively!" );
    }

    // // set up arrays of boundary conditions
    // std::vector<double> dummy;
    // // check if the user has given condition
    // if ( typeid(D2x2s_left) == typeid(dummy) ){
    //     // check if given list has correct size
    //     if ( d2x2s_left.size() != n ){
    //         throw std::length_error( "Left edge boundary array must have the same length as array of y coordinates!" );
    //     }
    //     // set internal variable
    //     //else{d2x2s_left = std::move(D2x2s_left);}
    // }
    // // otherwise, ignore and the boundary list will stay as array of zeros

    // // perform similar checks on all other boundary conditions
    // if ( typeid(D2x2s_right) == typeid(dummy) ){
    //     if ( d2x2s_right.size() != n ){
    //         throw std::length_error( "Right edge boundary array must have the same length as array of y coordinates!" );
    //     }
    //     //else{d2x2s_right = std::move(D2x2s_right);}
    // }

    // if ( typeid(D2y2s_bottom) == typeid(dummy) ){
    //     if ( d2y2s_bottom.size() != m ){
    //         throw std::length_error( "Bottom edge boundary array must have the same length as array of x coordinates!" );
    //     }
    //     //else{d2y2s_bottom = std::move(D2y2s_bottom);}
    // }

    // if ( typeid(D2y2s_top) == typeid(dummy) ){
    //     if ( d2y2s_top.size() != m ){
    //         throw std::length_error( "Top edge boundary array must have the same length as array of x coordinates!" );
    //     }
    //     // else{d2y2s_top = std::move(D2y2s_top);}
    // }

    // if ( typeid(D4x2y2s_corners) == typeid(dummy) ){
    //     if ( d4x2y2s_corners.size() != 4 ){
    //         throw std::length_error( "Corner boundary array must have length 4!" );
    //     }
    //     //else{d4x2y2s_corners = std::move(D4x2y2s_corners);}
    // }

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //! put this in a private method
    // initialise matrices
    // matrix of spline coefficients
    std::vector<std::vector<std::vector<double>>> S(n-1, std::vector<std::vector<double>>(m-1, std::vector<double>(16, 0.0)));
    // matrix of x derivatives at knots
    std::vector<std::vector<double>> zs_x(n, std::vector<double>(m, 0));

    // matrix of y derivatives at knots
    std::vector<std::vector<double>> zs_y(n, std::vector<double> (m, 0));

    // matrix of xy derivatives at knots
    std::vector<std::vector<double>> zs_xy(n, std::vector<double> (m, 0));

    // boundary vectors of xyy derivatives at top and bottom knots
    std::vector<double> zs_xyy_bottom(m, 0);
    std::vector<double> zs_xyy_top(m, 0);
    
    if ( verbose ){
        std::cout << "Constructing 1D splines along each grid line..." << std::endl;
    }

    // construct 1D splines along x and y to find first 
    // and cross derivatives at spline knots

    // construct splines along boundary rows
    CubicSpline spline_bottom(xs, zs[0], d4x2y2s_corners[0], d4x2y2s_corners[1]);
    CubicSpline spline_top(xs, zs[n-1], d4x2y2s_corners[3], d4x2y2s_corners[2]);

    std::cout<<"1"<<std::endl;
    // get spline coefficients
    // Change these to getters
    std::vector<std::vector<double>> S_bottom = spline_bottom.GetS();
    std::vector<std::vector<double>> S_top = spline_top.GetS();

    std::cout<<"2"<<std::endl;
    // calculate and store analytical x derivatives at boundary knots
    // the rest of the terms of the derivative are zero at the knots
    // due to (x - t) terms (see CubicSpline code)
    std::cout<<"3"<<std::endl;
    for (  int i = 0; i<zs_xyy_bottom.size()-1; i++){
        zs_xyy_bottom[i] = S_bottom[i][1];
    }
    zs_xyy_bottom[zs_xyy_bottom.size()-1] = S_bottom[S_bottom.size()-1][1] + 2*S_bottom[S_bottom.size()-1][2]*(xs[n-1] - xs[n-2]) + 3* S_bottom[S_bottom.size()-1][3]*std::pow((xs[n-1] - xs[n-2]),2);
    std::cout<<"4"<<std::endl;
    for (  int i = 0; i<zs_xyy_top.size()-1; i++ ){
        zs_xyy_top[i] = S_top[i][1];
    }
    zs_xyy_top[zs_xyy_top.size()-1] = S_top[S_top.size()-1][1] + 2* S_top[S_top.size()-1][2]*(xs[n-1] - xs[n-2]) + 3*S_top[S_top.size()-1][3]*std::pow((xs[n-1] - xs[n-2]),2);
    std::cout<<"5"<<std::endl;
    // loop over data rows
    for (  int i = 0; i<n; i++ ){
        
        // construct spline along ith row using the 
        // user input data (excluding the boundary rows)
        CubicSpline spline_x(xs, zs[i], d2x2s_left[i], d2x2s_right[i]);
        // get spline coefficients
        std::vector<std::vector<double>> S_x = spline_x.GetS();

        // calculate and store analytical x derivatives 
        // at knots along the ith row
        std::cout<<"1"<<std::endl; 
        for (  int j = 0; j<zs_x[i].size()-1; j++ ){
            zs_x[i][j] = S_x[j][1];
        }
        std::cout<<"1"<<std::endl; 
        zs_x[i][zs_x[i].size()-1] = S_x[S_x.size()-1][1] + 2*S_x[S_x.size()-1][2]*(xs[n-1] - xs[n-2]) + 3*S_x[S_x.size()-1][3]*std::pow((xs[n-1] - xs[n-2]),2);
    }
    std::cout<<"6"<<std::endl;
    // loop over data columns
    for (  int i = 0; i < m; i++ ){
        // construct spline along jth column using the calculated 
        // spline derivatives (excluding boundary columns)
        std::vector<double> temp; 
        for (  int j = 0; j<zs_x.size(); j++ ){
            temp.push_back(zs_x[j][i]);
        }
        CubicSpline spline_xy(ys, temp, zs_xyy_bottom[i], zs_xyy_top[i]);
        // get spline coefficients
        std::vector<std::vector<double>> S_xy = spline_xy.GetS();

        // calculate and store analytical xy derivatives
        // at knots along the jth column
        for (  int j = 0; j<zs_xy.size()-1; j++ ){
            zs_xy[j][i] = S_xy[j][1];
        }
        zs_xy[zs_xy.size()-1][i] = S_xy[S_xy.size()-1][1] + 2*S_xy[S_xy.size()-1][2]*(ys[m-1] - ys[m-2]) + 3*S_xy[S_xy.size()-1][3]*std::pow((ys[m-1] - ys[m-2]),2);
        
        // construct spline along jth column using the user input data
        // (excluding boundary columns since they have no further use)
        std::vector<double> temp2;
        for (  int j = 0; j<zs.size(); j++ ){
            temp2.push_back(zs[j][i]);
        }
        CubicSpline spline_y(ys, temp2, d2y2s_bottom[i], d2y2s_top[i]);

        // get spline coefficients
        std::vector<std::vector<double>> S_y = spline_y.GetS();

        // calculate and store analytical xy derivatives
        // at knots along the jth column
        for (  int j = 0; j<zs_y.size()-1; j++ ){
            zs_y[j][i] = S_y[j][1];
        }
        zs_y[zs_y.size()-1][i] = S_y[S_y.size()-1][1] + 2*S_y[S_y.size()-1][2]*(ys[m-1] - ys[m-2]) + 3*S_y[S_y.size()-1][3]*std::pow((ys[m-1] - ys[m-2]),2);
    }

    
    if ( verbose ){
        std::cout<<"Filling spline coefficient matrix..."<<std::endl;
    }

    // fill knot values and derivatives into matrix

    // loop over the outer dimensions of the matrix
    for (  int j=0; j<n-1; j++ ){
        for (  int i=0; i<m-1; i++ ){

            // array of knot values and knot derivatives
            // at each iteration (grid cell)
            std::vector<double>Z(16,0);

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

            if ( verbose ){
                if ( i == 0 && j == 0 ){
                    std::cout<<"Here is a matrix Z containing the function values and x, y, xy derivatives at the 4 knots of the first rectangle:"<<std::endl;
                    linalg::PrettyPrint(Z);
                }
            }

            // matrix A such that A * S[j][i] = Z
            std::vector<std::vector<double>> A(16, std::vector<double>(16,0));

            // loop over dimensions of A
            for (  int l=0; l<4; l++ ){
                for (  int k=0; k<4; k++ ){
                    // function value terms
                    //! Template recursion -- more efficient for integer powers
                    A[0][4*l+k] = std::pow(xs[i],k) * std::pow(ys[j],l);
                    A[1][4*l+k] = std::pow(xs[i+1],k) * std::pow(ys[j],l);
                    A[2][4*l+k] = std::pow(xs[i],k) * std::pow(ys[j+1],l);
                    A[3][4*l+k] = std::pow(xs[i+1],k) * std::pow(ys[j+1],l);

                    //initialise remaining terms to zero
                    A[4][4*l+k] = A[5][4*l+k] = A[6][4*l+k] = A[7][4*l+k] = 0;
                    A[8][4*l+k] = A[9][4*l+k] = A[10][4*l+k] = A[11][4*l+k] = 0;
                    A[12][4*l+k] = A[13][4*l+k] = A[14][4*l+k] = A[15][4*l+k] = 0;

                    if ( k > 0 ){
                        // x derivative terms
                        A[4][4*l+k] = k * std::pow(xs[i],k-1) * std::pow(ys[j],l);
                        A[5][4*l+k] = k * std::pow(xs[i+1],k-1) * std::pow(ys[j],l);
                        A[6][4*l+k] = k * std::pow(xs[i],k-1) * std::pow(ys[j+1],l);
                        A[7][4*l+k] = k * std::pow(xs[i+1],k-1) * std::pow(ys[j+1],l);

                        if ( l > 0 ){
                            // y derivative terms
                            A[8][4*l+k] = l * std::pow(xs[i],k) * std::pow(ys[j],l-1);
                            A[9][4*l+k] = l * std::pow(xs[i+1],k) * std::pow(ys[j],l-1);
                            A[10][4*l+k] = l * std::pow(xs[i],k) * std::pow(ys[j+1],l-1);
                            A[11][4*l+k] = l * std::pow(xs[i+1],k) * std::pow(ys[j+1],l-1);

                            // xy derivative terms
                            A[12][4*l+k] = k * l * std::pow(xs[i],k-1) * std::pow(ys[j],l-1);
                            A[13][4*l+k] = k * l * std::pow(xs[i+1],k-1) * std::pow(ys[j],l-1);
                            A[14][4*l+k] = k * l * std::pow(xs[i],k-1) * std::pow(ys[j+1],l-1);
                            A[15][4*l+k] = k * l * std::pow(xs[i+1],k-1) * std::pow(ys[j+1],l-1);
                        }
                    } else if ( l > 0 ){
                        // remaining y derivative terms
                        A[8][4*l+k] = l * std::pow(xs[i],k) * std::pow(ys[j],l-1);
                        A[9][4*l+k] = l * std::pow(xs[i+1],k) * std::pow(ys[j],l-1);
                        A[10][4*l+k] = l * std::pow(xs[i],k) * std::pow(ys[j+1],l-1);
                        A[11][4*l+k] = l * std::pow(xs[i+1],k) * std::pow(ys[j+1],l-1);
                    }
                }
            }
            
            if ( verbose ){
                if ( i == 0 && j == 0 ){
                    std::cout<<"Here is the corresponding matrix A:"<<std::endl;
                    linalg::PrettyPrint(A);
                }
            }
            
            // find determinant of A
            //! Delete? waste of computation when its not used
            double A_det = linalg::Determinant(A);

            if ( verbose ){
                if ( i == 0 && j == 0 ){
                    std::cout << "and its determinant:" << std::endl;
                    std::cout << A_det << std::endl;
                }
            }   

            // invert A
            //! possibly get eigen library to handle c++ linalg stuff 
            std::vector<std::vector<double>> A_inv = linalg::Inverse(A);

            if ( verbose ){
                if ( i == 0 and j == 0 ){
                    std::cout<<"and its inverse:"<<std::endl;
                    linalg::PrettyPrint(A_inv);
                }
            }

            // calculate the coefficient array for this iteration (grid cell)
            // meaning: S[j][i] = A^-1 * Z
            // S[j][i] = linalg::MatMul(A_inv, Z);
            S[j][i] = linalg::GuassElim(A, Z);

            if ( verbose ){
                if ( i == 0 and j == 0 ){
                    std::cout<<"The corresponding spline coefficients then read:"<<std::endl;
                    linalg::PrettyPrint(S[j][i]);
                }
            }
        }
    }
    // define callable attribute
    SetS(S);

    if ( verbose ){
        std::cout<<"Finally, this is the full spline coefficient matrix:"<<std::endl;
        //linalg::PrettyPrint(S);
        // find a nice way to print out 3D vector
    }
    
    // automatic print statement
    std::cout<<"Spline coefficients have been successfully calculated and stored!"<<std::endl;
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
    
    bool verbose = true;

    if ( verbose ){std::cout<<"in evaluateSpline():"<<std::endl;}

    double x = X;
    double y = Y;


    // initialise array storing function values
    double splines;

    // find indices before which each user input value should be
    // slotted into the arrays of knot coordinates such that the
    // arrays remain sorted; these indices (j, i) correspond to
    // the spline indices (i, j) at which to evaluate the spline
    // for each input point

    if ( verbose ){std::cout<< "Determining grid locations..." << std::endl;}

    // relevant indices
    // np.searchsorted uses binary search
    // side='right': xs[i] <= x_eval < xs[i+1] as required
    // indices given are those of the knot values before which each
    // input value should be slotted in so need to be shifted back by one
    int ixs,iys;
    ixs = std::upper_bound(xs.begin(), xs.end(), x) - xs.begin();
    iys = std::upper_bound(ys.begin(), ys.end(), y) - ys.begin();
    // if final input value equal to (or greater than) the final knot 
    // value, rightmost value evaluates to final index of knot array + 1 
    // so it needs to be shifted back by one more
    if (x >= xs[xs.size()-1] ){ixs -= 1;}
    
    if (y >= ys[ys.size()-1] ){iys -= 1;}

    if ( verbose ){
        std::cout<<"x index"<<std::endl;
        std::cout<<ixs;
        std::cout<<std::endl;
        std::cout<<"y index"<<std::endl;
        std::cout<<iys;
        std::cout<<std::endl;
        std::cout<<"Evaluating bicubic spline at each input coordinate..."<<std::endl;
    }

                
    // evaluate spline at (j, i)-th coordinates and store in array
    splines = (           _S[iys][iys][0]                + _S[iys][iys][1]  * x
                        + _S[iys][iys][2]  * std::pow(x,2)        + _S[iys][iys][3]  * std::pow(x,3)
                        + _S[iys][iys][4]         * y    + _S[iys][iys][5]  * x    * y
                        + _S[iys][iys][6]  * std::pow(x,2) * y    + _S[iys][iys][7]  * std::pow(x,3) * y
                        + _S[iys][iys][8]         * std::pow(y,2) + _S[iys][iys][9]  * x    * std::pow(y,2)
                        + _S[iys][iys][10] * std::pow(x,2) * std::pow(y,2) + _S[iys][iys][11] * std::pow(x,3) * std::pow(y,2)
                        + _S[iys][iys][12]        * std::pow(y,3) + _S[iys][iys][13] * x    * std::pow(y,3)
                        + _S[iys][iys][14] * std::pow(x,2) * std::pow(y,3) + _S[iys][iys][15] * std::pow(x,3) * std::pow(y,3) );
    return splines;
}

std::vector<std::vector<double>> BicubicSpline::evaluateSpline( std::vector<double> X, std::vector<double> Y ){
    /*
    Returns the value(s) of the interpolated function at x, y.

    params:  x:  x-coordinates(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
                y:  y-coordinate(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
    */
    
    bool verbose = false;

    if ( verbose ){std::cout<<"in evaluateSpline():"<<std::endl;}

    int m,n;
    std::vector<double> x = X;
    std::vector<double> y = Y;
    m = y.size();
    n = x.size();

    // initialise array storing function values
    std::vector<std::vector<double>> splines(m, std::vector<double>(n));

    // find indices before which each user input value should be
    // slotted into the arrays of knot coordinates such that the
    // arrays remain sorted; these indices (j, i) correspond to
    // the spline indices (i, j) at which to evaluate the spline
    // for each input point

    if ( verbose ){std::cout<< "Determining grid locations..." << std::endl;}

    // relevant indices
    // np.searchsorted uses binary search
    // side='right': xs[i] <= x_eval < xs[i+1] as required
    // indices given are those of the knot values before which each
    // input value should be slotted in so need to be shifted back by one
    std::vector<double> ixs(n);
    std::vector<double> iys(m);
    for ( int i=0; i<n; i++ ){
        ixs[i] = std::upper_bound(xs.begin(), xs.end(), x[i]) - xs.begin();
        iys[i] = std::upper_bound(ys.begin(), ys.end(), x[i]) - ys.begin();
    }
    // if final input value equal to (or greater than) the final knot 
    // value, rightmost value evaluates to final index of knot array + 1 
    // so it needs to be shifted back by one more
    if (x[x.size()-1] >= xs[xs.size()-1] ){ixs[ixs.size()-1] -= 1;}
        
    if (y[y.size()-1] >= ys[ys.size()-1] ){iys[iys.size()-1] -= 1;}

    if ( verbose ){
        std::cout<<"x indices"<<std::endl;
        for (const auto& val : ixs ){std::cout<<val<<" ";}
        std::cout<<std::endl;
        std::cout<<"y indices"<<std::endl;
        for (const auto& val : iys ){std::cout<<val<<" ";}
        std::cout<<std::endl;
        std::cout<<"Evaluating bicubic spline at each input coordinate..."<<std::endl;
    }

    // loop over input coordinates
    for ( int i=0; i<m; i++ ){
        for ( int j=0; j<n; j++ ){
                
            // evaluate spline at (j, i)-th coordinates and store in array
            splines[i][j] = (   _S[iys[i]][ixs[j]][0]                + _S[iys[i]][ixs[j]][1]  * x[j]
                                + _S[iys[i]][ixs[j]][2]  * std::pow(x[j],2)        + _S[iys[i]][ixs[j]][3]  * std::pow(x[j],3)
                                + _S[iys[i]][ixs[j]][4]         * y[i]    + _S[iys[i]][ixs[j]][5]  * x[j]    * y[i]
                                + _S[iys[i]][ixs[j]][6]  * std::pow(x[j],2) * y[i]    + _S[iys[i]][ixs[j]][7]  * std::pow(x[j],3) * y[i]
                                + _S[iys[i]][ixs[j]][8]         * std::pow(y[i],2) + _S[iys[i]][ixs[j]][9]  * x[j]    * std::pow(y[i],2)
                                + _S[iys[i]][ixs[j]][10] * std::pow(x[j],2) * std::pow(y[i],2) + _S[iys[i]][ixs[j]][11] * std::pow(x[j],3) * std::pow(y[i],2)
                                + _S[iys[i]][ixs[j]][12]        * std::pow(y[i],3) + _S[iys[i]][ixs[j]][13] * x[j]    * std::pow(y[i],3)
                                + _S[iys[i]][ixs[j]][14] * std::pow(x[j],2) * std::pow(y[i],3) + _S[iys[i]][ixs[j]][15] * std::pow(x[j],3) * std::pow(y[i],3) );
        }
    }
    if ( verbose ){
        std::cout<<"Here is the full array of function values:"<<std::endl;
        linalg::PrettyPrint(splines);
    }
        
    return splines;
}
// Setters
void BicubicSpline::SetS(std::vector<std::vector<std::vector<double>>> S){_S = S;}

// Getters
std::vector<std::vector<std::vector<double>>> BicubicSpline::GetS(){return _S;}
