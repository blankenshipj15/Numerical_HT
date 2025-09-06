/* Code is used to calculate specific values of the diffusion equation with a source
and the following conditions.
G.E. kT''=-50x^3, -T'=0 at x=0, T(x=6) = 100, domain is 6 units long
Created by Jesse Blankenship on 8/26/2025
Last edit: 2000 8/26/2025
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

int main(){
    // define Gamma as constant
    const double gamma = 2.0;
    int numRows;
    std::cout << "Enter the number of solution points: ";
    std::cin >> numRows;

    // define delta x for the domain
    double deltaX = 6.0 / static_cast<double>(numRows - 1);

    std::vector<double> x(numRows);
    std::vector<double> soln(numRows);

    // opens .csv file to write data for post-processing
    std::ofstream outputFile;
    outputFile.open("exactSoln.csv", std::ios::out); // this opens .csv file and will overwrite existing data

    // calculate solution at desired x values given the analytical solution
    for (int i = 0; i < numRows; ++i) {
        x[i] = static_cast<double>(i) * deltaX;
        soln[i] = (-50.0 * std::pow(x[i], 5)) / (gamma*20.0) + 9820.0;
        outputFile << static_cast<double>(x[i]) << "," << static_cast<double>(soln[i]) << "\n";

    }

    outputFile.close();

    return 0;
}
