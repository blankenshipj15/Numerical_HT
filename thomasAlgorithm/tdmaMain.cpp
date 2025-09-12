/* TDMA (Tridiagonal Matrix Algorithm, also known as Thomas Algorithm) main is the main file used for testing TDMA functions.
Tests cases are performed as desired to verify accuracy of functions.
Author: Jesse Blankenship
Last Updated: 9/12/2025
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "tdmaSolver.hpp"

int main(){

    // Define problem parameters
    const double gamma = 400; // conductivity
    const double tBcLeft = 500; // Dirichelet left BC, temperature
    const double tBcRight = 300; // Dirichelet right BC, temperature
    const double domainLen = 6; // total length of 1D domain
    const double nPoints = 3; // number of cell centroid temperatures calculated
    const double deltaX = domainLen/(nPoints-1);

    // Initialize the coefficient vectors for TDMA
    std::vector<double> a(nPoints), b(nPoints), c(nPoints), d(nPoints);
    std::vector<double> T(nPoints);

    a[0] = gamma/deltaX + 2*gamma/deltaX;
    a[1] = gamma/deltaX + gamma/deltaX;
    a[2] = 2*gamma/deltaX + gamma/deltaX;

    b[0] = gamma/deltaX;
    b[1] = gamma/deltaX;
    b[2] = 0;

    c[0] = 0;
    c[1] = gamma/deltaX;
    c[2] = gamma/deltaX;

    d[0] = 2*gamma*tBcLeft/deltaX;
    d[1] = 0;
    d[2] = 2*gamma*tBcRight/deltaX;

    // apply TDMA solver to the set of vectors
    T = tdmaSolver(a,b,c,d);

    for(int i = 0; i < nPoints; ++i){
        if(i == 0){
            std::cout << static_cast<double>(tBcLeft) << ", " << i*deltaX/2 << "\n";
        }
        std::cout << static_cast<double>(T[i]) << ", " << (i+1)*deltaX/2 << "\n";
        if(i == nPoints-1){
            std::cout << static_cast<double>(tBcRight) << ", " << (i+2)*deltaX/2 << "\n";
        }
    }

    return 0;
}