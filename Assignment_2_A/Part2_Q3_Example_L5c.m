close all; clear all; clc

%% motor parameters
VoltagePeakMin = 0.5;
VoltagePeakMax = 4;
N_RunningAverage =100;

%% Outrunner
J = 796.66972e-9; %from CAD: ABS - assume plastic is 42.9g
N_poles = 14; 
Data = csvread('RunDown_Outrunner_v1.csv',1,0);

%% extract run down data
Data = Data(10:end, :);
time = Data(:,1);
voltage = Data(:,2); 

%% filter data
B = 1/N_RunningAverage*ones(N_RunningAverage,1); 
voltage_filtered = filter(B,1,voltage);  

%% determine times for zero crossing and peaks
index_zc = find(voltage_filtered.*circshift(voltage_filtered,[1 1])<=0);
time_zc = time(index_zc);
index_peak = index_zc(1:end-1) + round(gradient(index_zc(1:end-1))/2);
time_peak = time(index_peak);
voltage_peak = voltage_filtered(index_peak);

%% trim points from before run down (start at 10th point after drops below threshold)
index_start = find(abs(voltage_peak)<VoltagePeakMax);
time_zc = time_zc(index_start(10):end);
time_peak = time_peak(index_start(10):end);
voltage_peak = voltage_peak(index_start(10):end);

%% trim points at end of curve to avoid spurious zero crossings
time_zc = time_zc(abs(voltage_peak)>VoltagePeakMin);
time_peak = time_peak(abs(voltage_peak)>VoltagePeakMin);
voltage_peak = voltage_peak(abs(voltage_peak)>VoltagePeakMin);

%% determine velocity
period = gradient(time_zc)*2; % two zerocrossings per full wave period
Hz_electrical = 1./period;
Hz_mechanical = Hz_electrical/(N_poles/2);
AngularVelocity = Hz_mechanical*2*pi;
VoltageRMS = abs(voltage_peak)/sqrt(2)/(sqrt(3)); % also convert from line-line to phase voltage
% VoltageRMS = abs(voltage_peak)/sqrt(2); % also convert from line-line to phase voltage (same for delta)

%% fit curve to determine torque constant
p_Ke = polyfit(AngularVelocity,VoltageRMS,1)
VoltageRMS_linear = polyval(p_Ke, AngularVelocity);
Kt = p_Ke(1)*3

%% polynominal fit for velocity and acceleration
TimePoly = linspace(min(time_zc)-2, max(time_zc)+0.49);  % take 2 seconds off the start time to simulate if the motor started faster
p_vel = polyfit(time_zc, AngularVelocity, 2);  % fit polynomial to AngularVelocity Data
p_acc = polyder(p_vel); % determine the acceleration as the gradient of the AngularVelocity Data

%% evaluate polynomial
VelocityPoly = polyval(p_vel, TimePoly);  % evaluate the polynomal for AngularVelocity
AccelerationPoly = polyval(p_acc, TimePoly);  % evaluate the polynomal for Acceleration

%% T = J*omega
T_spinning_loss = -J*AccelerationPoly;  % as per eq 53 in reference
P_spinning_loss = T_spinning_loss.*VelocityPoly;
Rpm = VelocityPoly * 60/(2*pi);  % define the rpm for plotting

Polynomial_SpinningLoss_Rpm = polyfit(Rpm, P_spinning_loss, 2);
Rpm_Plot = 0:600;
P_SpinningLoss_Plot = polyval(Polynomial_SpinningLoss_Rpm, Rpm_Plot);

%% plotting
figure
plot(time,voltage,time,voltage_filtered)
hold on
plot(time_zc,0,'x', time_peak, voltage_peak, 'x', 'markersize',10, 'linewidth', 2)
hold off
title('Run Down Voltage')
ylabel('Line to Line Voltage (V)')
xlabel('Time (sec)')
% 
figure
plot(AngularVelocity(2:end-1),VoltageRMS(2:end-1),'.',AngularVelocity(2:end-1), VoltageRMS_linear(2:end-1), 'linewidth', 2)
title('Voltage constant (Ke)')
ylabel('Phase Voltage (vrms)')
xlabel('Angular Velocity (rad/sec)')
axis([0 inf 0 inf])

figure
plot(Rpm, P_spinning_loss, Rpm_Plot, P_SpinningLoss_Plot)
title('spinning loss power')
xlabel('speed (rpm)')
ylabel('spinning loss (W)')
