#pragma once
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace linalg{
    template<typename T>
    void luDecomposition(std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& L, std::vector<std::vector<T>>& U, const int& n) {
        // Initialize L and U
        L = std::vector<std::vector<T>>(n, std::vector<T>(n, 0));
        U = std::vector<std::vector<T>>(n, std::vector<T>(n, 0));
        std::vector<std::vector<T>> temp = A;
        
        for (int i = 0; i < n; ++i) {
            L[i][i] = 1; // Set diagonal of L to 1
        }

        // Perform LU decomposition with partial pivoting
        for (int k = 0; k < n; ++k) {
            // Partial pivoting
            T maxVal = abs(A[k][k]);
            int maxRow = k;
            for (int i = k + 1; i < n; ++i) {
                if (std::abs(A[i][k]) > maxVal) {
                    maxVal = std::abs(A[i][k]);
                    maxRow = i;
                }
            }

            // Swap rows in A
            if (maxRow != k) {
                std::swap(A[maxRow], A[k]);
            }

            // Calculate U and L
            for (int j = k; j < n; ++j) {
                U[k][j] = A[k][j];
            }

            for (int i = k + 1; i < n; ++i) {
                L[i][k] = A[i][k] / U[k][k];
                for (int j = k; j < n; ++j) {
                    A[i][j] -= L[i][k] * U[k][j];
                }
            }
        }
        A = temp; //reset A as it points back to the value in main
    }
    
    template<typename T>
    std::vector<T> forwardSubstitution(const std::vector<std::vector<T>>& L, const std::vector<T>& b) {
        int n = L.size();
        std::vector<T> y(n);

        for (int i = 0; i < n; ++i) {
            y[i] = b[i];
            for (int j = 0; j < i; ++j) {
                y[i] -= L[i][j] * y[j];
            }
        }
        return y;
    }
    
    template<typename T>
    std::vector<T> backwardSubstitution(const std::vector<std::vector<T>>& U, const std::vector<T>& y) {
        int n = U.size();
        std::vector<T> x(n);

        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= U[i][j] * x[j];
            }
            x[i] /= U[i][i];
        }
        return x;
    }
    
    template<typename T>
    std::vector<T> luSolver(const std::vector<std::vector<T>>& L, const std::vector<std::vector<T>>& U, const std::vector<T>& b) {
        std::vector<T> y = forwardSubstitution(L, b);
        return backwardSubstitution(U, y);
    }
    
    template<typename T>
    std::vector<T> GuassElim(std::vector<std::vector<T>>& A, const std::vector<T>& b){
        //check square
        int N = A.size();
        if ( N!=A[0].size() ){
            throw std::length_error( "Error side M != side N, the matrix input must be square." );
        }
        

        // Initialise temp matrices
        std::vector<std::vector<T>> L,U;
        std::vector<T> sol;

        // Decompose then solve
        luDecomposition(A, L, U, N);
        sol = luSolver(L, U, b); 

        return sol;
    }
    
    template<typename T>
    T Determinant(std::vector<std::vector<T>>& A) {
        int N = A.size();
        std::vector<std::vector<T>> L, U;

        luDecomposition(A, L, U, N);

        T det = 1.0;
        for (int i = 0; i < N; ++i) {
            det *= U[i][i];
        }

        return det;
    }
    
    template<typename T>
    std::vector<std::vector<T>> Inverse(std::vector<std::vector<T>>& A) {
        int N = A.size();
        std::vector<std::vector<T>> L, U;

        luDecomposition(A, L, U, N);
        std::vector<std::vector<T>> inv(N, std::vector<T>(N));

        for (int i = 0; i < N; ++i) {
            std::vector<T> e(N, 0);
            e[i] = 1; // Create unit vector for column i
            std::vector<T> y = forwardSubstitution(L, e);
            std::vector<T> x = backwardSubstitution(U, y);
            for (int j = 0; j < N; ++j) {
                inv[j][i] = x[j];
            }
        }
        return inv;
    }
    
    template<typename T>
    void PrettyPrint(std::vector<std::vector<T>> A) {
        if ( A.empty() || A[0].empty() ) {
            throw std::length_error( "Error matrix is empty." );
        }
        // Find the maximum width needed for formatting
        int maxWidth = 0;
        for ( const auto& row : A ) {
            for ( const auto& num : row ) {
                maxWidth = std::max(maxWidth, static_cast<int>(std::to_string(num).size()));
            }
        }

        // Print the matrix with formatting
        for ( const auto& row : A ) {
            for ( const auto& num : row ) {
                std::cout << std::setw(maxWidth + 1) << num;  // +1 for spacing
            }
            std::cout << std::endl;
        }
    }
    
    template<typename T>
    void PrettyPrint(std::vector<T> b){
        if ( b.empty() ){
            throw std::length_error( "Error vector is empty." );
        }
        std::cout << "[ ";
        for (const auto& elem : b){
            std::cout << elem << " ";   
        }
        std::cout << "]" << std::endl;
    }
}