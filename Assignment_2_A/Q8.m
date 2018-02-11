%It has been calculated that the conductive area is Ac = 1.39mm2

%J in A/mm2
close all;
clear all;
clc;
J = [5:5:25];

I = J.*1.39;
Kt = 0.01519;%hand calculation and FEMM result
T = Kt.*I;

plot(J, T);
title('Torque vs Current Density');
xlabel('Current density J(A/mm^2)');
ylabel('Torque(Nm)');