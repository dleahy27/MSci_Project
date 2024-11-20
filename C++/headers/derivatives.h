#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>

struct Derivatives{
    std::vector<double> d2x2s_left; 
    std::vector<double> d2x2s_right;
    std::vector<double> d2y2s_bottom;
    std::vector<double> d2y2s_top;
    std::vector<double> d4x2y2s_corners;
};