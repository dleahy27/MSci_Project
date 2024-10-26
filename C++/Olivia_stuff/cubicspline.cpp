#include "cubicspline.h"

CubicSpline::CubicSpline(const std::vector<double>& Ts, const std::vector<double>& Ys, const double& D2t2_start = 0.0, const double& D2t2_end = 0.0) : N(Ts.size()), _S(N - 1, std::vector<double>(4, 0)), ts(Ts), ys(Ys), d2t2_start(D2t2_start), d2t2_end(D2t2_end)
{
    /*
    Finds cubic spline coefficients for given dataset.
    Each spline segment requires 4 coefficients.

    params:   ts: array of spline knots
                ys: array of function values at spline knots

    kwargs:   d2t2_start: second derivative at first knot
                            (default 0)
                d2t2_end:   second derivative at final knot
                            (default 0)
    */

    // define internal variables

    // check if the two lists are the same size
    if ( ys.size() != N ){
        throw std::length_error( "Arrays must have the same length!" );
    }

    // useful values for matrix calculation
    std::vector<double> hs(N-1,0);
    std::vector<double> bs(N-1,0);
    std::vector<double> us(N-1,0);
    std::vector<double> vs(N-1,0);

    // S''(ti) values
    std::vector<double> zs(N,0);

    // coefficients of each sub-spline,
    // where sub-splines Si have the form:
    // Si(x) = Ai + Bi(x-ti) + Ci(x-ti)^2 + Di(x-ti)^3
    // which is a Taylor expansion of Si about ti
    std::vector<std::vector<double>> S(N-1,std::vector<double>(4, 0));

    // some initial array elements
    zs[0] = d2t2_start;
    zs[N-1] = d2t2_end;

    // loop over size of arrays
    for (  int i = 0; i<N-1; i++ ){

        // define useful values for matrix calculation
        hs[i] = ts[i+1] - ts[i];
        bs[i] = (ys[i+1] - ys[i]) / hs[i];

        // forward elimination
        if ( i == 1 ){
            us[i] = 2 * (hs[i] + hs[i-1]);
            vs[i] = 6 * (bs[i] - bs[i-1]);
        } else {
            us[i] = 2 * (hs[i] + hs[i-1]) - ( std::pow(hs[i-1], 2) / us[i-1] );
            vs[i] = 6 * (bs[i] - bs[i-1]) - ( (hs[i-1] * vs[i-1]) / us[i-1] );
        }
    }
    // backwards loop over size of arrays
    for (  int i = N-2; i>-1; i -= 1 ){
        // back substitution
        if ( i == 0 ){
            S[i][0] = ys[i];
            S[i][1] = - (hs[i]*zs[i+1])/6 - (hs[i]*zs[i])/3 + bs[i];
            S[i][2] = zs[i]/2;
            S[i][3] = (zs[i+1] - zs[i])/(6*hs[i]);
        // find remaining coefficients
        } else{
        // the following works because zs[N-1] has been set to a constant
            zs[i] = ( vs[i] - (hs[i]*zs[i+1]) )/us[i];
            // find coefficients for each sub-spline
            S[i][0] = ys[i];
            S[i][1] = - (hs[i]*zs[i+1])/6 - (hs[i]*zs[i])/3 + bs[i];
            S[i][2] = zs[i]/2;
            S[i][3] = (zs[i+1] - zs[i])/(6*hs[i]);
        }
    }
    // define callable coefficient variable
    CubicSpline::SetS(S);
}

double CubicSpline::evaluateSpline(double x,std::vector<std::vector<double>>S, std::vector<double>ts, bool verbose=false){
    /*
    Returns the value(s) of the interpolated function at x.

    params:  x:  value(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
    */
    
    // initialise variables (define empty lists because np.empty()
    // needs to know the dtype in advance, but dtype could be
    // list of bools or bool for conditions and function for
    // subsplines)
    int N = GetN();
    std::vector<bool> conditions; 
    std::vector<std::function<double(double,  int)>> subsplines;
    bool cond;

    // loop over spline coefficients
    for (  int i = 0; i<N-1; i++ ){            
        // define condition for ith sub-spline: each condition is a 
        // boolean list defining the location along the spline of 
        // each user input value
        // at end of spline, value can also be equal to final knot
        if ( i == N-2 ){
            cond = (x >= ts[i] && x <= ts[i+1]);
            conditions.push_back(cond);
        }
        // otherwise
        else{
            cond = (x >= ts[i] && x < ts[i+1]);
            conditions.push_back(cond);
        }
        // define ith sub-spline
        auto subspline = [S, ts](double y,  int j){return ( S[j][0] + S[j][1]*(y - ts[j]) 
                                        + S[j][2]*std::pow(y - ts[j], 2) 
                                        + S[j][3]*std::pow(y - ts[j], 3) );};
        subsplines.push_back(subspline);
    }
        
    // piecewise define the full spline and evaluate 
    // at input values
    return piecewise(x, conditions, subsplines);
}

std::vector<double> CubicSpline::evaluateSpline(std::vector<double> x,std::vector<std::vector<double>>S, std::vector<double>ts, bool verbose=false){
    /*
    Returns the value(s) of the interpolated function at x.

    params:  x:  value(s) at which to evaluate function.
                    User can provide either a single input 
                    value or an array of inputs.
    */
    
    // initialise variables (define empty lists because np.empty()
    // needs to know the dtype in advance, but dtype could be
    // list of bools or bool for conditions and function for
    // subsplines)
    int N = GetN();
    std::vector<std::vector<bool>> conditions; 
    std::vector<std::function<double(double,  int)>> subsplines;

    // loop over spline coefficients
    for (  int i = 0; i<N-1; i++ ){           
        // define condition for ith sub-spline: each condition is a 
        // boolean list defining the location along the spline of 
        // each user input value

        std::vector<bool> cond;
        //! can think of optimising below
        for (  int j = 0; j<x.size(); j++ ){  
            // at end of spline, value can also be equal to final knot
            if ( i == N-2 ){
                cond.push_back(x[j] >= ts[i] && x[j] <= ts[i+1]);
            }
            // otherwise
            else{
                cond.push_back(x[j] >= ts[i] && x[j] < ts[i+1]);
            }
            
        }
        conditions.push_back(cond);
        // define ith sub-spline
        auto subspline = [S, ts](double y,  int j){return ( S[j][0] + S[j][1]*(y - ts[j]) 
                                        + S[j][2]*std::pow(y - ts[j], 2) 
                                        + S[j][3]*std::pow(y - ts[j], 3) );};
        subsplines.push_back(subspline);
    }           
    // piecewise define the full spline and evaluate 
    // at input values
    return piecewise(x, conditions, subsplines);
}

// Functions
double CubicSpline::piecewise(double x, std::vector<bool> conds, std::vector<std::function<double(double,  int)>> lambda){
    double val;
    for (  int i=0; i<(conds.size()-1); i++ ){
        if (conds[i]){
            val = lambda[i](x, i);
            return val; // for only one value no need to continue for loop
        }
    }
    return 0.0;
    // return error

}

std::vector<double> CubicSpline::piecewise(std::vector<double> x, std::vector<std::vector<bool>> conds, std::vector<std::function<double(double,  int)>> lambda){
    std::vector<double> vals;
    for (  int i=0; i<conds.size()-1; i++ ){
        //! same optimisation here maybe better way than looping twice
        for (  int j=0; j<x.size(); j++ ){
            if (conds[i][j]){
                vals.push_back(lambda[i](x[j], i));
            }
        }
    }
    // have if statement here for error
    if ( vals.size() != x.size() ){
        //error
        //return smth
    }
    return vals;
}

// Getters
std::vector<std::vector<double>> CubicSpline::GetS(){return _S;}
int CubicSpline::GetN(){return N;}

// Setters
void CubicSpline::SetS(std::vector<std::vector<double>> S){_S = S;}

