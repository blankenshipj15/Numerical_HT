/* TDMA (Tridiagonal Matrix Algorithm, also known as Thomas Algorithm) main is the main file used for testing TDMA functions.
Tests cases are performed as desired to verify accuracy of functions.
Author: Jesse Blankenship
Last Updated: 9/13/2025
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "tdmaSolver.hpp"

int main(){

    // Define problem parameters
    // *var*1 represents first case tested (1D diffusion with constant temperature BCs)
    // *var*2 is second case (1D diffusion with constant teperature and flux BCs)
    const double gamma = 400; // conductivity
    const double tBcLeft1 = 500; // Dirichlet left BC, temperature
    const double tBcRight1 = 300; // Dirichlet right BC, temperature
    const double tBcLeft2 = 500; // Dirichlet left BC, temperature
    const double fluxBcRight2 = 3000; // Neumann right BC, flux
    const double domainLen = 6; // total length of 1D domain
    const double nPoints = 3; // number of cell centroid temperatures calculated
    const double deltaX = 2; // distance between cell centroids

    // Initialize the coefficient vectors for TDMA with test case 1
    std::vector<double> a1(nPoints), b1(nPoints), c1(nPoints), d1(nPoints);
    std::vector<double> T1(nPoints);

    a1[0] = gamma/deltaX + 2*gamma/deltaX;
    a1[1] = gamma/deltaX + gamma/deltaX;
    a1[2] = 2*gamma/deltaX + gamma/deltaX;

    b1[0] = gamma/deltaX;
    b1[1] = gamma/deltaX;
    b1[2] = 0;

    c1[0] = 0;
    c1[1] = gamma/deltaX;
    c1[2] = gamma/deltaX;

    d1[0] = 2*gamma*tBcLeft1/deltaX;
    d1[1] = 0;
    d1[2] = 2*gamma*tBcRight1/deltaX;

    // Initialize the coefficient vectors for TDMA with test case 2
    std::vector<double> a2(nPoints), b2(nPoints), c2(nPoints), d2(nPoints);
    std::vector<double> T2(nPoints);

    a2[0] = gamma/deltaX + 2*gamma/deltaX;
    a2[1] = gamma/deltaX + gamma/deltaX;
    a2[2] = gamma/deltaX;

    b2[0] = gamma/deltaX;
    b2[1] = gamma/deltaX;
    b2[2] = 0;

    c2[0] = 0;
    c2[1] = gamma/deltaX;
    c2[2] = gamma/deltaX;

    d2[0] = 2*gamma*tBcLeft1/deltaX;
    d2[1] = 0;
    d2[2] = -fluxBcRight2;

    // apply TDMA solver to the set of vectors and return the solution vector of temperatures
    T1 = tdmaSolver(a1,b1,c1,d1);
    std::cout << "Case 1: \n";
    for(int i = 0; i < nPoints; ++i){
        if(i == 0){
            std::cout << static_cast<double>(tBcLeft1) << "\n";
        }
        std::cout << static_cast<double>(T1[i]) << "\n";
        if(i == nPoints-1){
            std::cout << static_cast<double>(tBcRight1) << "\n" << "\n";
        }
    }

    T2 = tdmaSolver(a2,b2,c2,d2);
    std::cout << "Case 2: \n";
    for(int i = 0; i < nPoints; ++i){
        if(i == 0){
            std::cout << static_cast<double>(tBcLeft1) << "\n";
        }
        std::cout << static_cast<double>(T2[i]) << "\n";
        if(i == nPoints-1){
            std::cout << static_cast<double>((((2*gamma/deltaX)+(gamma/deltaX))*
                                                T2[nPoints-1]-(gamma/deltaX)*T2[nPoints-2])/(2*gamma/deltaX))
                                                 << "\n" << "\n";
        }
    }

    return 0;
}