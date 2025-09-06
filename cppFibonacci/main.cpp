/*
main.cpp created to execute a Fibonacci sequence generator.
Created by: Jesse Blankenship
Date: 08-20-2025
*/

#include <iostream>
#include"fib.hpp"

int main() {
    int i = 0;
    while(i == 0){
        std::cout << "Enter the desired fibonacci number index (it should be an integer): ";
        int n;
        std::cin >> n;
        int fibNumber = fibonacci(n);
        std::cout << "Fibonacci number at index " << n << " is: " << fibNumber << std::endl;
        std::cout << "Would you like to calculate another Fibonacci number? (y/n) ";
        char response;
        std::cin >> response;
        if(response != 'y'){
            ++i;
        }
        std::cout << std::endl << std::endl;
    }
    return 0;
}

