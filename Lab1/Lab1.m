clc;
clear;
close all;
addpath 'include'

%% Parameters
alfa = 1.4;
beta = 1.8;
gamma = 1.05;
m = 13.5;
J = 0.37;

%% Matrices 
% Define the system matrices

% Q1. Considering the state vector and the input:
%   x = [x y theta vx vy omega]^T
%   u = [Fx Fy Nz]

% Q2 We can define x. = Ax + Bu
A = [ % cinematics part:
      0 0 0 1 0 0; 
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      % dynamics part:
      0 0 0 -alfa/m 0 0;
      0 0 0 0 -beta/m 0;
      0 0 0 0 0 -gamma/J ]; 

B = [ 0 0 0;
      0 0 0;
      0 0 0;
      1/m 0 0;
      0 1/m 0;
      0 0 1/J ];

%% Simulation
dt = 0.01; % integration step
t = 0:dt:10; % time vector

%% Input
u = zeros(3, length(t)); % input
u(1,:) = 3 * ones(1,length(t));
u(2,:) = 3 * ones(1,length(t));
u(3,:) = 2 * ones(1,length(t));

%% Initial state
x0 = [1;2;3;0;0;0]; % initial state

%% Simulation using analytical formula
% Q3. Estimate x(t) thru analytical formula
x = SimulateSysAnalytical(A, B, u, t, x0, dt);

% position x1 and x2 increase linearly in time -> coherent with model

%% Simulation using Euler integration
% Q4. Estimate x(t) thru euler integration
xEuler = SimulateSysEuler(A, B, u, t, x0, dt);

% same result as Q3: position x1 and x2 increase linearly in time

%% Plot
% Q5. Plot  state estimate

% Plot the state evolution (analytical)
figure;
plot(t, x(1,:), 'r', 'DisplayName', 'x_1');
hold on;
plot(t, x(2,:), 'b', 'DisplayName', 'x_2');
% plot other states...
xlabel('t[s]');
ylabel('x');
title('System state (analytical formula)');
legend;
grid on;

% Plot the state evolution (Euler)
figure;
plot(t, xEuler(1,:), 'r', 'DisplayName', 'x_1');
hold on;
plot(t, xEuler(2,:), 'b', 'DisplayName', 'x_2');
% plot other states...
xlabel('t[s]');
ylabel('x');
title('System state (Euler integration)');
legend;
grid on;

% Plot difference in percentage of Euler - analytical
figure;
plot(t, (xEuler(1,:) - x(1,:))./x(1,:) * 100, 'r', 'DisplayName', 'x_1');
hold on;
plot(t, (xEuler(2,:) - x(2,:))./x(2,:) * 100, 'b', 'DisplayName', 'x_2');
% plot other states...
xlabel('t[s]');
ylabel('Difference percentage %');
title('State estimation difference percentage % (Euler - analytical)');
legend;
grid on;

%% Q6. Difference between the euler - analytical estimation

% - At the beginning, the percentage error is zero, which indicates that the two simulations (analytical and Euler) coincide perfectly at the initial time
% - % Subsequently, the error becomes slightly negative, reaching a minimum value of approximately −0.2% for x1 and −0.3% for x2

% This means that 
% - the Euler method slightly underestimates the analytical solution (explicit integration methods)
% - the error remains very small and stable, showing that, with an integration step of d=0.01
% the Euler method provides a good approximation of the analytical solution, while avoiding the calculation of the matrix exponential and the convolution integral.


%% Q7. Possible Output

% Considering LTI system: 
%   x.=Ax+Bu 
%   y=Cx 
%   with state x = [x y theta vx vy omega]^T

% The output is a linear combination of this variables, that depends on the sensors we have to measure this varibles. 
