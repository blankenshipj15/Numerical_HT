clear;
clc;

fileID = fopen('exactSoln.csv','r','n');
solnData = cell2mat(textscan(fileID,'%f%f', 'Delimiter', ','));
fclose(fileID);

plot(solnData(:,1), solnData(:,2))