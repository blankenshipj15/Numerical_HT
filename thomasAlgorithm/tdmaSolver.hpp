/* A function that takes in vector coefficients for linear equations arranged in a tridiagonal matrix
form and returns the solution vector using the TDMA method.
Author: Jesse Blankenship
Last Updated: 9/12/2025
*/

#include <iostream>
#include <vector>
std::vector<double> tdmaSolver(
    const std::vector<double> &a, // the a coefficients of the main diagonal (phi at i)
    const std::vector<double> &b, // the b coefficients for phi at i+1
    const std::vector<double> &c, // the c coefficients for phi at i-1
    const std::vector<double> &d  // all costants for the equation defining phi at i
){
    // define the P, Q, and phi (solution) vectors
    const int nPoints = a.size();
    std::vector<double> P(nPoints), Q(nPoints), phi(nPoints);

    P[0] = b[0]/a[0];
    Q[0] = d[0]/a[0];
    for(int i = 1; i < nPoints; ++i){
        P[i] = b[i]/(a[i]-c[i]*P[i-1]);
        Q[i] = (d[i] + c[i]*Q[i-1])/(a[i] - c[i]*P[i-1]);
    }
    
    phi[nPoints-1] = Q[nPoints-1];
    for(int i = nPoints-2; i >= 0; --i){
        phi[i] = P[i]*phi[i+1] + Q[i];
    }

    return phi;
};