/* tdma2Dmain is the main file used for testing the TDMA function in a 2D heat conduction problem.
Boundary conditions are all Dirichlet, conductivity is assumed constant, the domain is rectangular.
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
    const double tBcLeft = 293; // Dirichlet left BC, temperature
    const double tBcRight = 373; // Dirichlet right BC, temperature
    const double tBcTop = 373; // Dirichlet left BC, temperature
    const double tBcBottom = 293; // Neumann right BC, flux
    const double width = 1; // total length of 1D domain
    const double height = 1; // number of cell centroid temperatures calculated
    const double nXpoints = 3;
    const double nYpoints = 3;
    const double deltaX = width/(nXpoints-1);
    const double deltaY = height/(nYpoints-1);

    // Initialize the coefficient matrices for TDMA on 2D domain
    std::vector<double> a(nXpoints, nYpoints), b(nXpoints, nYpoints), c(nXpoints, nYpoints), d(nXpoints, nYpoints);
    std::vector<double> T(nXpoints, nYpoints);

    return 0;
}