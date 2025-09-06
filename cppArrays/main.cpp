
#include <iostream>
#include <vector>
#include <cmath>


int main(){

    std::vector<double> myArray(3);

    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            myArray[i]= 1;
        } 
    }
    
    return 0;
}