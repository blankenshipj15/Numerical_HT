clc;
clear;
% conductivity value
gamma = 400;

% numerical result for dirichlet-dirichlet
% results from tdmaMain.cpp
T_num1 = [500; 
        466.667;
        400;
        333.333;
        300];
% numerical result for dirichlet-neumann
% results from tdmaMain.cpp
T_num2 = [500;
        492.5;
        477.5;
        462.5;
        455];

x = [0;1;3;5;6];
T_analytic1 = -200*gamma.*x./(6*gamma) + 500;
T_analytic2 = -3000.*x./gamma + 500;

figure; 
subplot(1,2,1)
plot(x, T_num1, "-O")
hold on
plot(x,T_analytic1, "--")
legend(["Numerical-TDMA", "Analytic"])
title("Dirichlet-Dirichlet")
hold off

subplot(1,2,2)
plot(x, T_num2, "-O")
hold on
plot(x,T_analytic2, "--")
legend(["Numerical-TDMA", "Analytic"])
title("Dirichlet-Neumann")
hold off

