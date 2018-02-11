close all; clear all; clc

%% motor parameters
VoltagePeakMin = 0.5;
VoltagePeakMax = 4;
N_RunningAverage =100;

%% Outrunner
N_poles = 14; 
Data = csvread('Copy_of_2018-01-29_170841.csv',1,0);

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
VoltageRMS = abs(voltage_peak)/sqrt(2)/(sqrt(3)); % also convert from line-line to phase voltage (for star)
% VoltageRMS = abs(voltage_peak)/sqrt(2); % also convert from line-line to phase voltage (same for delta)

%% fit curve to determine torque constant
p_Ke = polyfit(AngularVelocity,VoltageRMS,1)
VoltageRMS_linear = polyval(p_Ke, AngularVelocity);
Kt = p_Ke(1)*3

%% plotting
figure
plot(time,voltage,time,voltage_filtered)
hold on
plot(time_zc,0,'x', time_peak, voltage_peak, 'x', 'markersize',10, 'linewidth', 2)
hold off
title('Run Down Voltage')
ylabel('Line to Line Voltage (V)')
xlabel('Time (sec)')

figure
plot(AngularVelocity(2:end-1),VoltageRMS(2:end-1),'.',AngularVelocity(2:end-1), VoltageRMS_linear(2:end-1), 'linewidth', 2)
title('Voltage constant (Ke)')
ylabel('Phase Voltage (vrms)')
xlabel('Angular Velocity (rad/sec)')
axis([0 inf 0 inf])

