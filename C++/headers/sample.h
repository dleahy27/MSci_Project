#pragma once 
#include <cmath>
#include <vector>
#include <numbers>
#define PI 3.14159265359

namespace sample{
    template <typename T>
    std::vector<T> linspace(T min, T max, int n){
        std::vector<T> x(n);
        // fill grid knot arrays i.e. linspace
        for (  int i = 0; i<n; i++){
            if ( i == n-1 ){
                x[i] = max;
            } else {
                T val = min + i*(max - min)/(n-1); 
                x[i] = val;
            }  
        }
        return x;
    }

    template <typename T>
    std::vector<T> logspace(T min, T max, int n, float base = 10.0){
        std::vector<T> logx(n);
        std::vector<T> x(n);
        // fill grid knot arrays i.e. linspace
        x = linspace(min,max,n);
        for (  int i = 0; i<n; i++){
            logx[i] = std::pow(base,x[i]); 
        }
        return logx;
    }

    template <typename T>
    std::vector<T> chebyshev(T min, T max, int n, unsigned int kind = 1){
        std::vector<T> x;
        x.reserve(n);

        // affine transformation parameters
        double affine_x0 = (max + min)/2; // offset
        double affine_a = (max - min)/2; // factor

        // fill grid knot arrays
        if ( kind == 1 ){
            for ( int i = 0; i < n; i++ ){
                x[i] = std::cos(PI*(2*i + 1)/(2*n));// chebyshev nodes of the first kind in interval (min,max)
            }
        } else if ( kind == 2 ){
            for ( int i = 0; i < n; i++ ){
                x[i] = std::cos(PI*i/(n - 1));//chebyshev nodes of the second kind in interval [min,max]
            }
        } else{
            // raise error
            throw std::invalid_argument( "Only first kind and second kind Chebyshev nodes exist." );
        }

        return x;
    }
}