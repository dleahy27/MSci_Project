#include "../headers/cubicspline.h"
#include <algorithm>

// fast square function
inline double square(double x){return x*x;}

// fast cube function
inline double cube(double x){return x*x*x;}

CubicSpline::CubicSpline(const std::vector<double>& Ts, const std::vector<double>& Ys, double D2t2_start = 0.0, double D2t2_end = 0.0) : N(Ts.size()), _S(N - 1, std::vector<double>(4)), ts(Ts), ys(Ys), d2t2_start(D2t2_start), d2t2_end(D2t2_end)
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
    std::vector<double> hs(N-1);
    std::vector<double> bs(N-1);
    std::vector<double> us(N-1);
    std::vector<double> vs(N-1);

    // S''(ti) values
    std::vector<double> zs(N);

    // coefficients of each sub-spline,
    // where sub-splines Si have the form:
    // Si(x) = Ai + Bi(x-ti) + Ci(x-ti)^2 + Di(x-ti)^3
    // which is a Taylor expansion of Si about ti
    std::vector<std::vector<double>> S(N-1,std::vector<double>(4));

    // some initial array elements
    us[0] = 0;
    vs[0] = 0;
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
        } 
        if ( i > 1 ){
            us[i] = 2 * (hs[i] + hs[i-1]) - (( hs[i-1]*hs[i-1]) / us[i-1] );
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
    SetS(S);
}

double CubicSpline::evaluateSpline(double t, const std::vector<std::vector<double>>& S, const std::vector<double>& ts){
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

    int it; 
    std::function<double(double,  int)> subspline;

    it = std::upper_bound(ts.begin(), ts.end(), t) - ts.begin() - 1;

    if (t >= ts[ts.size()-1] ){it -= 1;}

    // define ith sub-spline
    subspline = [S, ts](double y,  int j){return ( S[j][0] + S[j][1]*(y - ts[j]) 
                                        + S[j][2]*square(y - ts[j]) 
                                        + S[j][3]*cube(y - ts[j]) );};
        
    // piecewise define the full spline and evaluate 
    // at input values
    return piecewise(t, it, subspline);
}

std::vector<double> CubicSpline::evaluateSpline(const std::vector<double>& t, const std::vector<std::vector<double>>& S, const std::vector<double>& ts){
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

    std::function<double(double,  int)> subspline;

    std::vector<int> its;

    // loop over spline coefficients
    for (  int i = 0; i<t.size(); i++ ){           
        // define condition for ith sub-spline: each condition is a 
        // boolean list defining the location along the spline of 
        // each user input value
        its[i] = std::upper_bound(ts.begin(), ts.end(), t[i]) - ts.begin() - 1;

        if (t[t.size()-1] >= ts[ts.size()-1] ){its[its.size()-1] -= 1;}

        
    }

    // define sub-spline
    subspline = [S, ts](double y,  int j){return ( S[j][0] + S[j][1]*(y - ts[j]) 
                                    + S[j][2]*square(y - ts[j]) 
                                    + S[j][3]*cube(y - ts[j]) );};           
    // piecewise define the full spline and evaluate 
    // at input values
    return piecewise(t, its, subspline);
}

// Functions
double CubicSpline::piecewise(const double t, const int it, std::function<double(double,  int)>& lambda){
    return lambda(t, it);
    
}

std::vector<double> CubicSpline::piecewise(const std::vector<double>& t, const std::vector<int>& its, std::function<double(double,  int)>& lambda){
    std::vector<double> vals;
    for (  int i=0; i<its.size(); i++ ){
        vals[i] = lambda(t[i], its[i]);
        }
    return vals;
}

// Getters
std::vector<std::vector<double>> CubicSpline::GetS(){return _S;}

// Setters
void CubicSpline::SetS(const std::vector<std::vector<double>>& S){_S = S;}

