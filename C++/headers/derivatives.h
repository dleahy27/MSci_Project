#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>

/// @brief Stores the corner and boundary derivatives of the grid.
///
/// This struct stores various partial derivatives that are needed for the bicubic algorithm:
/// the second partial derivatives of the grid with respect to y at the top edge and bottom edge of the grid,
/// the second partial derivatives of the grid with respect to x at the left edge and right edge of the grid,
/// the fourth partial derivative with respect to x and y (2x2y) at the 4 corner points of the grid.
struct Derivatives{
    std::vector<double> d2x2s_left; 
    std::vector<double> d2x2s_right;
    std::vector<double> d2y2s_bottom;
    std::vector<double> d2y2s_top;
    std::vector<double> d4x2y2s_corners;
};