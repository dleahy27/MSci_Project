#include "../headers/funcs.h"

inline constexpr double power(double base, int exponent) {
    if (exponent == 0) {
        return 1;
    } else if (exponent > 0) {
        return base * power(base, exponent - 1);
    } else {
        return 1 / power(base, -exponent);
    }
}

// random generator
std::vector<double> funcs::random_uniform(double start, double end, unsigned int N){
    std::vector<double> sol;
    sol.reserve(N);

    std::default_random_engine gen;
    std::uniform_real_distribution<double> distribution(start,end);

    for (int i = 0; i<N; i++){
        sol.push_back(distribution(gen));
    }
    return sol;
}
    
double funcs::gauss(const double& x, const double& y){return ( std::cos(y)*std::sin(x) + 1.5 )*std::exp(-(x*x + y*y)/6);}

std::vector<std::vector<double>> funcs::gauss(const std::vector<double>& x, const std::vector<double>& y){
    int n,m;

    m = x.size();
    n = y.size();
    
    std::vector<std::vector<double>> z(n, std::vector<double> (m));
    
    for (  int i=0; i<n; i++ ){
        for (  int j=0; j<m; j++  ){
            z[i][j] = ( std::cos(y[i]) * std::sin(x[j]) + 1.5 )*std::exp(-(x[j]*x[j] + y[i]*y[i])/6);
        }
    }
    
    return z;
}

std::vector<std::vector<double>> funcs::trig(const std::vector<double>& x, const std::vector<double>& y){
    int n,m;

    m = x.size();
    n = y.size();
    
    std::vector<std::vector<double>> z(n, std::vector<double> (m));
    
    for (  int i=0; i<n; i++ ){
        for (  int j=0; j<m; j++  ){
            z[i][j] = ( std::cos(y[i]) * std::sin(x[j]) + 1.5 );
        }
    }
    
    return z;
}

std::vector<std::vector<double>> funcs::pdfLike(const std::vector<double>& x, const std::vector<double>& y){
    int n,m;

    m = x.size();
    n = y.size();
    
    std::vector<std::vector<double>> z(n, std::vector<double> (m));
    
    for (  int i=0; i<n; i++ ){
        for (  int j=0; j<m; j++  ){
            z[i][j] = ( (std::pow(x[j], -3) + std::pow(1-x[j],5) -1)* std::log10(y[i])) ;
        }
    }
    
    return z;
}

std::vector<double> funcs::d2dx2_gauss(const double x, const std::vector<double>& y){
    std::vector<double> sol(y.size());
    double cos_x = std::cos(x);
    double sin_x = std::sin(x);

    for ( int i=0; i<y.size(); i++ ){
        double exp_xy = std::exp((1.0/6.0)*(-x*x - y[i]*y[i]));
        double cos_y = std::cos(y[i]);
        
        sol[i] = -(2.0/3.0)*exp_xy*x*cos_x*cos_y 
                - exp_xy*cos_y*sin_x 
                + (-(1.0/3.0)*exp_xy + (1.0/9.0)*exp_xy*x*x)*(1.5 + cos_y*sin_x);
    }
    return sol;
}

std::vector<double> funcs::d2dy2_gauss(const std::vector<double>& x, const double y){
    std::vector<double> sol(x.size());
    double cos_y = std::cos(y);
    double sin_y = std::sin(y);

    for ( int i=0; i<x.size(); i++ ){
        double exp_xy = std::exp((1.0/6.0)*(-x[i]*x[i] - y*y));
        double sin_x = std::sin(x[i]);

        sol[i] = -exp_xy*cos_y*sin_x 
                + (-(1.0/3.0)*exp_xy + (1.0/9.0)*exp_xy*y*y)*(1.5 + cos_y*sin_x) 
                + (2.0/3.0)*exp_xy*y*sin_x*sin_y;
        }

    return sol;
}

std::vector<double> funcs::d4dx2dy2_gauss(const std::vector<double>& X, const std::vector<double>& Y){
    std::vector<double> sol(4);
    std::vector<double> temp(4);

    // calculate d4dx2dy2 in loop
    for ( int i=0; i<2; i++ ){
        for ( int j=0; j<2; j++ ){
        
            double x = X[i];
            double y = Y[j];

            double exp_xy = std::exp((1.0/6.0)*(-x*x - y*y));
            double cos_x = std::cos(x);
            double cos_y = std::cos(y);
            double sin_x = std::sin(x);
            double sin_y = std::sin(y);

        temp[j + 2*i] = (2.0/3.0)*exp_xy*x*cos_x*cos_y 
                - (2.0/3.0)*x*(-(1.0/3.0)*exp_xy + (1.0/9.0)*exp_xy*y*y)*cos_x*cos_y 
                + exp_xy*cos_y*sin_x 
                - (-(1.0/3.0)*exp_xy + (1.0/9.0)*exp_xy*x*x)*cos_y*sin_x 
                + ((1.0/3.0)*exp_xy - (1.0/9.0)*exp_xy*y*y)*cos_y*sin_x 
                + (1.0/3.0*((1.0/3.0)*exp_xy - (1.0/9.0)*exp_xy*y*y) 
                + (1.0/9.0)*x*x*(-(1.0/3.0)*exp_xy + (1.0/9.0)*exp_xy*y*y))*(1.5 + cos_y*sin_x) 
                - (4.0/9.0)*exp_xy*x*y*cos_x*sin_y 
                - (2.0/3.0)*exp_xy*y*sin_x*sin_y 
                - 2.0*((1.0/9.0)*exp_xy*y - (1.0/27.0)*exp_xy*x*x*y)*sin_x*sin_y;
        }
    }

    // rearrange to bicubic format (counter-clockwise: [0,0] [1,0] [1,1] [0,1])
    sol[0] = temp[0];
    sol[3] = temp[1];
    sol[1] = temp[2];
    sol[2] = temp[3];

    return sol;
}

Derivatives funcs::deriv_gauss(const std::vector<double>& x, const std::vector<double>& y){
    Derivatives ds;

    int m = x.size()-1;
    int n = y.size()-1;

    ds.d2x2s_left = d2dx2_gauss(x[0], y);
    ds.d2x2s_right = d2dx2_gauss(x[m], y);
    ds.d2y2s_bottom = d2dy2_gauss(x, y[0]);
    ds.d2y2s_top = d2dy2_gauss(x, y[n]);
    ds.d4x2y2s_corners = d4dx2dy2_gauss(std::vector<double> {x[0],x[m]},std::vector<double> {y[0],y[n]});

    return ds;
}

std::vector<double> funcs::d1dx_trig(const double x, const std::vector<double>& y){
    std::vector<double> sol(y.size());

    double cosx = std::cos(x);  // Precompute cosx for efficiency

    for (int i = 0; i<y.size(); i++) {
        sol[i] = cosx * std::cos(y[i]); 
    }
    return sol;
}

std::vector<double> funcs::d1dy_trig(const std::vector<double>& x, const double y){
    std::vector<double> sol(x.size());

    double siny = std::sin(y);  // Precompute cosx for efficiency

    for (int i = 0; i<x.size(); i++) {
        sol[i] =  -siny * std::sin(x[i]); 
    }
    return sol;
}

std::vector<double> funcs::d2dx2_trig(const double x, const std::vector<double>& y){
    std::vector<double> sol(y.size());
    double sinx = std::sin(x);

    for ( int i=0; i<y.size(); i++ ){      
        sol[i] = -sinx*std::cos(y[i]);
    }
    return sol;
}

std::vector<double> funcs::d2dy2_trig(const std::vector<double>& x, const double y){
    std::vector<double> sol(x.size());
    double cosy = std::cos(y);
    
    for ( int i=0; i<x.size(); i++ ){      
        sol[i] = -std::sin(x[i])*cosy;
    }
    return sol;
}

std::vector<double> funcs::d4dx2dy2_trig(const std::vector<double>& x, const std::vector<double>& y){
    std::vector<double> sol(4);
    std::vector<double> temp(4);

    // calculate d4dx2dy2 in loop
    for ( int i=0; i<2; i++ ){
        for ( int j=0; j<2; j++ ){temp[j + 2*i] = std::cos(y[j]) + std::sin(x[i]);}
    }

    // rearrange to bicubic format (counter-clockwise: [0,0] [1,0] [1,1] [0,1])
    sol[0] = temp[0];
    sol[3] = temp[2];
    sol[1] = temp[3];
    sol[2] = temp[4];

    return sol;
}

Derivatives funcs::deriv_trig(const std::vector<double>& x, const std::vector<double>& y){
    Derivatives ds;

    int m = x.size()-1;
    int n = y.size()-1;

    ds.d2x2s_left = d2dx2_trig(x[0], y);
    ds.d2x2s_right = d2dx2_trig(x[m], y);
    ds.d2y2s_bottom = d2dy2_trig(x, y[0]);
    ds.d2y2s_top = d2dy2_trig(x, y[n]);
    ds.d4x2y2s_corners = d4dx2dy2_trig(std::vector<double> {x[0],x[m]},std::vector<double> {y[0],y[n]});

    return ds;
}

std::vector<double> funcs::d1dx_pdfLike(const double x, const std::vector<double>& y){
    std::vector<double> sol(y.size());

    double numerator_constant = -5*power(-1 + x, 4) - (3 / power(x, 4));
    double log_10 = std::log(10);  // Precompute log(10) for efficiency

    for (int i = 0; i<y.size(); i++) {
        double numerator = numerator_constant * std::log(y[i]); 
        sol[i] = (numerator / log_10); 
    }
    return sol;
}

std::vector<double> funcs::d1dy_pdfLike(const std::vector<double>& X, const double y){
    std::vector<double> sol(X.size());

    double denominator = y * std::log(10);

    for (int i = 0; i<X.size(); i++) {
        double numerator = -1 + power(1 - X[i], 5) + ( 1 / power(X[i], 3));
        sol[i] = (numerator / denominator);
    }
    return sol; 
}

std::vector<double> funcs::d2dx2_pdfLike(const double x, const std::vector<double>& y){
    std::vector<double> sol(y.size());

    double numerator_constant = -20 * power(-1 + x, 3) + 12 / power(x, 5);
    double log_10 = std::log(10);  // Precompute log(10) for efficiency

    for (int i = 0; i<y.size(); i++) {
        double numerator = numerator_constant * std::log(y[i]); 
        sol[i] = (numerator / log_10); 
    }
    return sol;
}

std::vector<double> funcs::d2dy2_pdfLike(const std::vector<double>& X, const double y){
    std::vector<double> sol(X.size());

    double denominator = power(y, 2) * std::log(10);

    for (int i = 0; i<X.size(); i++) {
        double numerator = -(-1 + power(1 - X[i], 5) + 1 / power(X[i], 3));
        sol[i] = (numerator / denominator);
    }
    return sol; 
}

std::vector<double> funcs::d4dx2dy2_pdfLike(const std::vector<double>& X, const std::vector<double>& Y){
    std::vector<double> sol(4);
    std::vector<double> temp(4);

    double log_10 = std::log(10);

    // calculate d4dx2dy2 in loop
    for ( int i=0; i<2; i++ ){
        for ( int j=0; j<2; j++ ){
        double numerator = 20 * power(-1 + X[i], 3) - 12 / power(X[i], 5);
        double denominator = power(Y[j], 2) * log_10;

        temp[j + 2*i] = (numerator/denominator);
        }
    }
    // rearrange to bicubic format (counter-clockwise: [0,0] [1,0] [1,1] [0,1])
    sol[0] = temp[0];
    sol[3] = temp[1];
    sol[1] = temp[2];
    sol[2] = temp[3];

    return sol;
}

Derivatives funcs::deriv_pdfLike(const std::vector<double>& x, const std::vector<double>& y){
    Derivatives ds;

    int m = x.size()-1;
    int n = y.size()-1;

    ds.d2x2s_left = d2dx2_pdfLike(x[0], y);
    ds.d2x2s_right = d2dx2_pdfLike(x[m], y);
    ds.d2y2s_bottom = d2dy2_pdfLike(x, y[0]);
    ds.d2y2s_top = d2dy2_pdfLike(x, y[n]);
    ds.d4x2y2s_corners = d4dx2dy2_pdfLike(std::vector<double> {x[0],x[m]},std::vector<double> {y[0],y[n]});

    return ds;
}

void funcs::outputTrigDerivatives( const std::vector<double>& x, const std::vector<double>& y, const std::string& filename ){
    const unsigned int m = x.size();
    const unsigned int n = y.size();

    std::vector<std::vector<double>> d1x_temp(m, std::vector<double>(n));
    std::vector<std::vector<double>> d2x_temp(m, std::vector<double>(n));
    std::vector<std::vector<double>> d1y_temp(n, std::vector<double>(m));
    std::vector<std::vector<double>> d2y_temp(n, std::vector<double>(m));
    
    std::ofstream myfile;
    myfile.open("../outputs/"+filename);
    myfile << "d1x,d1y,d2x,d2y"<<std::endl;

    for ( int i = 0; i<m; i++ ){
        d1x_temp[i] = d1dx_trig(x[i],y);
        d2x_temp[i] = d2dx2_trig(x[i],y);
    }

    for ( int j = 0; j<n; j++ ){
        d1y_temp[j] = d1dy_trig(x,y[j]);
        d2y_temp[j] = d2dy2_trig(x,y[j]);
        
    }

    for( int i = 0; i<n; i++ ){
        for( int j = 0; j<m; j++ ){
            myfile<<d1x_temp[j][i]<<","<<d1y_temp[i][j]<<","<<d2x_temp[j][i]<<","<<d2y_temp[i][j]<<std::endl;
        }
    }
    myfile.close();
}

void funcs::outputPdfDerivatives( const std::vector<double>& x, const std::vector<double>& y, const std::string& filename ){
    const unsigned int m = x.size();
    const unsigned int n = y.size();

    std::vector<std::vector<double>> d1x_temp(m, std::vector<double>(n));
    std::vector<std::vector<double>> d2x_temp(m, std::vector<double>(n));
    std::vector<std::vector<double>> d1y_temp(n, std::vector<double>(m));
    std::vector<std::vector<double>> d2y_temp(n, std::vector<double>(m));
    
    std::ofstream myfile;
    myfile.open("../outputs/"+filename);
    myfile << "d1x,d1y,d2x,d2y"<<std::endl;

    for ( int i = 0; i<m; i++ ){
        d1x_temp[i] = d1dx_pdfLike(x[i],y);
        d2x_temp[i] = d2dx2_pdfLike(x[i],y);
    }

    for ( int j = 0; j<n; j++ ){
        d1y_temp[j] = d1dy_pdfLike(x,y[j]);
        d2y_temp[j] = d2dy2_pdfLike(x,y[j]);
        
    }

    for( int i = 0; i<n; i++ ){
        for( int j = 0; j<m; j++ ){
            myfile<<d1x_temp[j][i]<<","<<d1y_temp[i][j]<<","<<d2x_temp[j][i]<<","<<d2y_temp[i][j]<<std::endl;
        }
    }
    myfile.close();
}