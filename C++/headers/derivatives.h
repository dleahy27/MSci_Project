#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>

/// @brief Stores the corner and boundary derivatives of the grid.
///
/// This struct stores various partial derivatives that are needed for the bicubic algorithm:
/// the second partial derivatives of the grid with respect to Q2 at the top edge and bottom edge of the grid,
/// the second partial derivatives of the grid with respect to x at the left edge and right edge of the grid,
/// the fourth partial derivative with respect to x and Q2 (2x2Q2) at the 4 corner points of the grid.
struct Derivatives{
    std::vector<double> d2xfd2xs_left; 
    std::vector<double> d2xfd2xs_right;
    std::vector<double> d2xfd2Q2s_bottom;
    std::vector<double> d2xfd2Q2s_top;
    std::vector<double> d4xfd2xd2Q2s_corners;
};