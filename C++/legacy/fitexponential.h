#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>

class ExponentialFitter {
public:
    // Constructor that initializes with separate x and y data vectors
    ExponentialFitter(const std::vector<double>& x_values, const std::vector<double>& y_values)
        : x_data(x_values), y_data(y_values) {
        if (x_data.size() != y_data.size()) {
            throw std::invalid_argument("x and y data vectors must be the same size.");
        }
    }

    // Method to fit the data and return the parameters [a, b, c]
    std::vector<double> fit(double learning_rate = 0.001, int max_iters = 10000, double tolerance = 1e-6) {
        double a = 1.0, b = 1.0, c = 1.0;  // Initial guesses for parameters

        for (int iter = 0; iter < max_iters; ++iter) {
            double grad_a = 0.0;
            double grad_b = 0.0;
            double grad_c = 0.0;

            // Compute gradients for a, b, and c
            for (size_t i = 0; i < x_data.size(); ++i) {
                double x = x_data[i];
                double y = y_data[i];
                double prediction = a * std::pow(x, b) + c;
                double error = prediction - y;

                grad_a += 2 * error * std::pow(x, b);         // ∂error/∂a
                grad_b += 2 * error * a * std::pow(x, b) * std::log(x);  // ∂error/∂b
                grad_c += 2 * error;                          // ∂error/∂c
            }

            // Update parameters
            a -= learning_rate * grad_a / x_data.size();
            b -= learning_rate * grad_b / x_data.size();
            c -= learning_rate * grad_c / x_data.size();

            // Check for convergence
            if (std::abs(grad_a) < tolerance && std::abs(grad_b) < tolerance && std::abs(grad_c) < tolerance) {
                break;
            }
        }

        return {a, b, c};
    }

    // Method to predict y for a given x based on fitted parameters
    double predict(double x, const std::vector<double>& params) const {
        double a = params[0];
        double b = params[1];
        double c = params[2];
        return a * std::pow(x, b) + c;
    }

private:
    std::vector<double> x_data;
    std::vector<double> y_data;
};