/*
Fibonacci sequence generator header file.
Created by: Jesse Blankenship
Date: 08-20-2025
*/

int fibonacci(int n) {
    if (n <= 1) {
        return n; // Base cases: fib(0) = 0, fib(1) = 1
    }

    int prev1 = 1; // Represents F(i-1)
    int prev2 = 0; // Represents F(i-2)
    int current;
    for (int i = 2; i <= n; ++i) {
        current = prev1 + prev2;
        prev2 = prev1;
        prev1 = current;
        // std::cout << current << "\n"; // used to debug
    }
    return current;
}