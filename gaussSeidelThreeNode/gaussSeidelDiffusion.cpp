/* This code is used to calculate the solution to the diffusion equation with a spatially 
varying source term in 1D. A zero flux condition is placed at the left boundary and
a constant value is set at the right boundary. Finite volume is used to discretize.
Gauss-Seidel is used to iteratively solve.
Author: Jesse Blankenship
Last Updated: 9/1/2025
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

int main(){

    // define constant gamma and delta x for problem givens
    const double gamma = 2.0;
    const double deltaX = 2.0;

    // initialize phi solution vectors and b vector for Ax=b
    std::vector<double> phiGuess = {{10000}, {8000}, {5000}};
    std::vector<double> phiNew(3);
    std::vector<double> bVector(3);

    // build b vector using integrated source term 50x^3
    bVector = {{-(50.0/4.0)*std::pow(2,4)}, {-(50.0/4.0)*(std::pow(4,4)-std::pow(2,4))},
               {-(50.0/4.0)*(std::pow(6,4)-std::pow(4,4))}};
    
    // define convergence criteria
    const double tolerance = 1.0e-5;
    const int maxIter = 1e3;
    int iter = 0;
    double errorMax;
    double error = 500.0;

    while(error >= tolerance){

        // leftmost cell centroid
        phiNew[0] = (bVector[0] - (gamma/deltaX) * phiGuess[1]) / (-gamma/deltaX);

        // since only three cell centroids considered, no loop is required for interior cell centroids
        phiNew[1] = (bVector[1] - (gamma/deltaX) * phiNew[0] - (gamma/deltaX) *phiGuess[2]) / (-2*gamma/deltaX);

        // rightmost cell centroid
        phiNew[2] = (bVector[2] - (gamma*2/deltaX) * 100 - (gamma/deltaX) * phiNew[1]) / ((-2*gamma/deltaX) - gamma/deltaX);

        // calculate error
        std::vector<double> errorVector(3);
        for(int i = 0; i < 3; ++i){
            errorVector[i] = std::abs(phiGuess[i]-phiNew[i]);
        }

        std::cout << "Iteration: " << iter <<", Phi 1: " << phiNew[0] << ", Phi 2: " << phiNew[1] << ", Phi 3: " << phiNew[2] << "\n";
        // find the maximum error in the error vector
        auto errorMaxPtr = std::max_element(errorVector.begin(), errorVector.end());
        errorMax = *errorMaxPtr;
        error = errorMax;
        
        // update solution vector with newly calculated values
        phiGuess = phiNew;
        iter = ++iter;

        if(iter > maxIter){
            break;
        }
    }

    // print solution data to screen
    std::cout << "The solution converged in " << iter << " iterations. \n";
    std::cout << "The solution vector is: \n";
    for(int i = 0; i < 3; ++i){
        std::cout << "Phi " << i+1 << " = " << static_cast<double>(phiNew[i]) << "\n";
    }

    return 0;
}