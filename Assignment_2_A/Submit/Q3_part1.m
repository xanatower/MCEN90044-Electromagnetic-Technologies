close all; clear all; clc

%% input data
NumberPoles = 14; 

%% extract FEMM data
Data = csvread('OutrunnerVideoExample_Part3.csv',1,0);
Angle = Data(:,1);
Torque = Data(:,2);
FluxLinkage = Data(:,3:5);
Current = Data(:,6:8);

%% calculations
K = circshift(FluxLinkage, [-round(length(Angle)*(90)/360),0])*NumberPoles/2;  % differentiate in frequency domain
T = K.*Current;

%% plotting
figure
plot(Angle, FluxLinkage)
hold on
plot(Angle,  K,'--')
hold off
legend('\lambda_a', '\lambda_b', '\lambda_c', 'k_a', 'k_b', 'k_c')
xlabel('electrical angle (degrees)')
ylabel('Wb & Vs')
axis([0 360 -inf inf])

figure
plot(Angle,  K,'--')
hold on
plot(Angle, Current)
hold off
legend('k_a', 'k_b', 'k_c', 'i_a', 'i_b', 'i_c')
xlabel('electrical angle (degrees)')
ylabel('Vs & A')
axis([0 360 -inf inf])


figure
plot(Angle,  T,'--')
hold on
plot(Angle, sum(T,2),Angle , Torque)
hold off
legend('T_a', 'T_b', 'T_c', 'T_{reconstructed}', 'T_{FEMM}')
xlabel('electrical angle (degrees)')
ylabel('Nm')
axis([0 360 -inf inf])

FluxLinkageRMS = max(max(FluxLinkage))/sqrt(2);
AverageTorque = mean(Torque);


Ke = FluxLinkageRMS*NumberPoles/2
Kt = Ke*3 %reconstructed, flux linkage from FEMM

%Irms
Irms = max(max(Current))/sqrt(2);
Kt_2 = AverageTorque/Irms 


