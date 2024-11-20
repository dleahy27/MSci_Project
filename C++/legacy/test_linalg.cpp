#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include "linalg.h"

int main(int argc, char* argv[]){
    // Initialise variables
    std::vector<std::vector<double>> mat1 = {{5,3,9},{-2,3,-1},{-1,-4,5}};
    std::vector<std::vector<double>> mat2 = {{1.0,2.0},{2.0,1.0}};
    std::vector<double> vec1 = {-1,-2,1};

    //std::vector<int> *Piv = new std::vector<int>;

    double an1 = 61.0/187.0;
    std::cout << an1 << std::endl;
    double an2 = -92.0/187.0;
    double an3 = -24.0/197.0;

    std::vector<double> an_solution = {an1,an2,an3};
    float anal_det = 187.0;

    std::vector<double> num_solution = linalg::GuassElim(mat1, vec1);
    double det = linalg::Determinant(mat1);
    std::vector<std::vector<double>> inv1 = linalg::Inverse(mat1);
    std::vector<std::vector<double>> inv2 = linalg::Inverse(mat2);

    std::cout << "The analytical solution is " << std::endl;
    linalg::PrettyPrint(an_solution);
    std::cout << std::endl << "The numerical solution is " << std::endl;
    linalg::PrettyPrint(num_solution);
    std::cout << std::endl;
    std::cout << "The numerical determinant is " << det << std::endl;
    std::cout << std::endl << "The analytical determinant is " << anal_det << std::endl;
    std::cout << std::endl << "The matrix is " << std::endl;
    linalg::PrettyPrint(mat1);
    std::cout << std::endl << "The numerical inverse is " << std::endl;
    linalg::PrettyPrint(inv1);
    std::cout << std::endl;
    std::cout << std::endl << "Another matrix is " << std::endl;
    linalg::PrettyPrint(mat2);
    std::cout << std::endl << "Another numerical inverse is " << std::endl;
    linalg::PrettyPrint(inv2);
    std::cout << std::endl << std::endl;


    return 0;
}