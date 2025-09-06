#include<iostream>
#include<Eigen/Eigen>

using namespace Eigen;


int main(){

    // Setting values of matrix with loop
    Matrix <double, 3, 3> matrixA;
    matrixA.setOnes();
    std::cout << matrixA << std::endl;
    std::cout << "\n";

    for(int i = 0; i<3; ++i){
        for(int j = 0; j<3; ++j){

            matrixA(i,j) = i*j;

        }
    }

    std::cout << matrixA << std::endl;
    std::cout << "\n";

    Matrix <double, 4, 4> matrixB;
    Matrix <double, 4, 1> matrixC;

    // Extracting diagonal values of one matrix and generating a column vector of those values
    matrixB = 4.0 * matrixB.setOnes();
    matrixC = matrixB.diagonal();
    std::cout << matrixC << std::endl;
    std::cout << "\n";


    // Identity matrix
    MatrixXd matrixD(10, 10);
    matrixD.setIdentity();
    std::cout << matrixD << std::endl;
    std::cout << "\n";
    std::cout << "Number of rows: " << matrixD.rows() << " Number of columns: " << matrixD.cols() << std::endl;
    std::cout << "\n";

    // Practice with linspaced function in Eigen
    VectorXd vectorA;
    vectorA = vectorA.setLinSpaced(100, 0, 1);
    std::cout << vectorA << std::endl;
    std::cout << "\n";


    return 0;
}