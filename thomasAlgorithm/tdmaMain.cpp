/* TDMA (Tridiagonal Matrix Algorithm, also known as Thomas Algorithm) main is the main file used for testing TDMA functions.
Tests cases are performed as desired to verify accuracy of functions.
Author: Jesse Blankenship
Last Updated: 9/10/2025
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

int main(){

    const double gamma = 400;
    const double tBcLeft = 500;
    const double tBcRight = 300;
    const double domainLen = 6;
    const double nPoints = 5;
    const double deltaX = domainLen/(nPoints-1);

    // Initialize the coefficient vectors for TDMA
    std::vector<double> a(nPoints), b(nPoints), c(nPoints), d(nPoints);
    std::vector<double> P(nPoints), Q(nPoints);
    std::vector<double> T(nPoints);

    b[0] = gamma/deltaX;
    b[1] = gamma/deltaX;
    b[2] = 0;

    c[0] = 0;
    c[1] = gamma/deltaX;
    c[2] = gamma/deltaX;

    for(int i = 0; i < nPoints; ++i){
        a[i] = b[i] + c[i];
    }

    d[0] = tBcLeft;
    d[1] = 0;
    d[2] = tBcRight;

    for(int i = 0; i < nPoints; ++i){
        if(i == 0){
            P[i] = b[i]/a[i];
            Q[i] = d[i]/a[i];
        } else{ 
            std::cout << (a[i]-c[i]*P[i-1]) << "\n";
            P[i] = b[i]/(a[i]-c[i]*P[i-1]);
            Q[i] = (d[i] + c[i]*Q[i-1])/(a[i] - c[i]*P[i-1]);
        }

    }
        
    return 0;
}