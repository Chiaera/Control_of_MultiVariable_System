clear;
clc;
close all;

%% set SCENARIO
% 1 = Q2. Luenberg observer, only position measured (x, y, yaw)
% 2 = Q3. Trivial observer,  full state measured

scenario = 2;

%% SET SYSTEM PARAMETERS
alfa = 1.4;
beta = 1.8;
gamma = 1.05;
m = 13.5;
J = 0.37;

% considering state x = [x y theta vx vy omega]^T
A = [ 0,0,0,1,0,0;  % state-space matrix
      0,0,0,0,1,0;
      0,0,0,0,0,1;
      0,0,0,-alfa/m,0,0;
      0,0,0,0,-beta/m,0;
      0,0,0,0,0,-gamma/J ];

xf = [-1;2;3;0;0;0]; % final state

B = [0 0 0; % input matrix
    0 0 0; 
    0 0 0;
    1/m 0 0;
    0 1/m 0;
    0 0 1/J];
usize = min(size(B));

T = 25; % simulation time
dt = 0.01; % step integration for Euler method

%% Q1. OBSERVABILITY

% correct scenario
switch scenario
    case 1
        disp('Q2: only position measured'); % Only position sensors: y = [x y yaw]^T
        C = zeros(3,6); 
        C(1:3,1:3) = diag([1,1,1]); 
    case 2
        disp('Q3: full state measured'); % All 6 states are directly measured
        C = eye(6); 
end

O = [ C;  % Kalman observability matrix
      C*A;
      C*A^2;
      C*A^3;
      C*A^4;
      C*A^5 ];

% check observability with rank(O)
rankO = rank(O);
disp(['rank(O) = ' num2str(rankO)]);

if rankO == size(A,1) % If rank(O) == 6: full-rank
    disp("Full-rank: system fully OBSERVABLE");
    
    % Q1 COMMENT:
    % Knowing 
    % - time evolution position x 
    % - physical responde of the model
    % ==> possible recontruct vx without a speed sensor
else
    disp("System is NOT observable");
end

%% Q2. LUENBERG OBSERVER
syms tau real

tic
eA = expm(A*(T-tau));
dG = eA*B*B'*eA';
G = vpa(int(dG,tau,0,T));

rankG = rank(G);
disp(['rank(G) = ' num2str(rankG)]);

if rankG == 6
    iG = inv(G);
else
    iG = pinv(G);
    warning("System has rank " + rankG + "! It is not fully controllable.")
end

u(tau) = B' * eA' * iG * xf; % IN to bring x_0 to x_f

toc
time = [0:dt:T];
time = time';
LL = length(time);

u_time = u(time);
uvalue = [eval(u_time{1}), eval(u_time{2}), eval(u_time{3})];

x0 = [0 0 0 0 0 0]'; % we start from the origin
x = zeros(6,LL); % this is only needed to log the results and plot them
xHat = zeros(6,LL); % this is only needed to log the results and plot them
x(:,1) = x0;

% SET OBSERVER PARAMETERS
% Set Ko (observability)
Tau = 1;
Ko = zeros(6,3);
xi = 1;
wn = 1;

% gain on position
kx1 = 2*xi*wn - alfa/m;
kx2 = 2*xi*wn - beta/m;
kx3 = 2*xi*wn - gamma/J;
% gain on velocity
kx4 = wn^2 - kx1*alfa/m;
kx5 = wn^2 - kx2*beta/m;
kx6 = wn^2 - kx3*gamma/J;

% define Ko [6x3]
Ko(1:3,1:3) = diag([kx1,kx2,kx3]);
Ko(4:6,1:3) = diag([kx4,kx5,kx6]);

%% SIMULATION
xHat(:,1) = x0 + [1, -1, 1.5, 0.5, -0.5, -1.5]'; % initial error on state extimation
y = zeros(min(size(C)),LL);
yHat = zeros(min(size(C)),LL);
yHat(:,1) = C * xHat(:,1);

for k = 1:LL-1
    % Luenberg equation: xHat. = A*xHat + B*u + K_o*(y - yHat)

    % Control
    u_k = uvalue(k,:)';
    xdot_k = A*x(:,k) + B*u_k;
    x(:,k+1) = x(:,k) + xdot_k * dt; % Euler integration

    % Measurement update
    y(:,k) = C * x(:,k);

    % Observer based on scenario
    switch scenario
        case 1  % Dynamic Luenberger Observer
            dx_est = A * xHat(:,k) + B * u_k + Ko * (y(:,k) - yHat(:,k));
            xHat(:,k+1) = xHat(:,k) + dx_est * dt; % Euler error on extimate state
            yHat(:,k+1) = C * xHat(:,k+1); % prediction

        case 2  % Trivial Observer: estimate coincides immediately with measurement
            xHat(:,k+1) = y(:,k); 
            yHat(:,k+1) = y(:,k);
    end
end

%% PLOTS
% Plot states
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Set figure to full screen
sgtitle("State") 
% Subplot 1: X Position
subplot(2,3,1);
plot(time, x(1,:),'linewidth',3);
hold on;
plot(time, ones(size(time)) * xf(1), 'linewidth',3);
xlabel('Time[s]');
ylabel('X[m]');
legend('State','Reference');
grid on;
title('X');

% Subplot 2: Y Position
subplot(2,3,2);
plot(time, x(2,:), 'linewidth',3);
hold on;
plot(time, ones(size(time)) * xf(2), 'linewidth',3);
xlabel('Time[s]');
ylabel('Y[m]');
legend('State','Reference');
grid on;
title('Y');

% Subplot 3: Theta
subplot(2,3,3);
plot(time, x(3,:), 'linewidth',3);
hold on;
plot(time, ones(size(time)) * xf(3), 'linewidth',3);
xlabel('Time[s]');
ylabel('Yaw[rad]');
legend('State','Reference');
grid on;
title('Yaw');

% Subplot 4: v_x
subplot(2,3,4);
plot(time, x(4,:), 'linewidth',3);
hold on;
plot(time, ones(size(time)) * xf(4), 'linewidth',3);
xlabel('Time[s]');
ylabel('v_x[m/s]');
legend('State','Reference');
grid on;
title('v_x');

% Subplot 5: v_y
subplot(2,3,5);
plot(time, x(5,:), 'linewidth',3);
hold on;
plot(time, ones(size(time)) * xf(5), 'linewidth',3);
xlabel('Time[s]');
ylabel('v_y[m/s]');
legend('State','Reference');
grid on;
title('v_y');

% Subplot 6: Yaw rate
subplot(2,3,6);
plot(time, x(6,:),'linewidth',3);
hold on;
plot(time, ones(size(time)) * xf(6), 'linewidth',3);
xlabel('Time[s]');
ylabel('Yaw rate[rad/s]');
legend('State','Reference');
grid on;
title('Yaw rate');

%% State stimation error
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Set figure to full screen
sgtitle("State estimation error") 
% Subplot 1: X Position
subplot(2,3,1);
plot(time, x(1,:) - xHat(1,:),'linewidth',3);
xlabel('Time[s]');
ylabel('X[m]');
grid on;
title('X');

% Subplot 2: Y Position
subplot(2,3,2);
plot(time, x(2,:) - xHat(2,:), 'linewidth',3);
xlabel('Time[s]');
ylabel('Y[m]');
grid on;
title('Y');

% Subplot 3: Theta
subplot(2,3,3);
plot(time, x(3,:) - xHat(3,:), 'linewidth',3);
xlabel('Time[s]');
ylabel('Yaw[rad]');
grid on;
title('Yaw');

% Subplot 4: v_x
subplot(2,3,4);
plot(time, x(4,:) - xHat(4,:), 'linewidth',3);
xlabel('Time[s]');
ylabel('v_x[m/s]');
grid on;
title('v_x');

% Subplot 5: v_y
subplot(2,3,5);
plot(time, x(5,:) - xHat(5,:), 'linewidth',3);
xlabel('Time[s]');
ylabel('v_y[m/s]');
grid on;
title('v_y');

% Subplot 6: Yaw rate
subplot(2,3,6);
plot(time, x(6,:) - xHat(6,:), 'linewidth',3);
xlabel('Time[s]');
ylabel('Yaw rate[rad/s]');
grid on;
title('Yaw rate');


figure;
plot(time,uvalue(:,1), 'linewidth',3)
hold on
plot(time,uvalue(:,2), 'linewidth',3)
plot(time,uvalue(:,3), 'linewidth',3)
legend('u_1','u_2','u_3')
xlabel('Time [s]')
ylabel('Inputs')
title('Input ')
grid on;

% Assuming you have time vector t corresponding to your data
figure;
scatter(x(1,:), x(2,:), [], time, 'filled'); % Plot with color mapping to time
colormap(jet); % Set colormap
cb = colorbar; % Add colorbar
cb.Label.String = 'time [s]'; % Add label to colorbar
xlabel('x [m]');
ylabel('y [m]');
title('Trajectory'); % Corrected title
grid on;

disp(['Final state error: ', num2str((x(:,end) - xf)')]);

disp(['Final est error: ', num2str((xHat(:,k+1) - x(:,k+1))')]);



%% PLOT COMMENTS

%% PLOT COMMENTS
% Q2 - Dynamic Luenberger Observer (Only Position Measured)
    
    %   State Estimation Error:
    %       It is possible to notice the transient difference between x and xHat.
    %       All the graphs start at the correct initial offset values considered in the text.
    %       In less than 10s, all profiles converge asymptotically to zero without oscillations.
    %       Indeed, with xi = 1 and wn = 1, the final estimation error is negligible (order of 10^-10), meaning it is numerically null.
    
    %   State:
    %       The state has the exact same profile as Lab 2, which confirms that the control law is working properly and independently of the observer dynamics.
   
    %   Trajectory:
    %       The trajectory does not show any uncertainty; the profile is a clean and straight diagonal line.
    
    %   Input:
    %       The input shows a continuous and smooth profile, confirming an optimal minimum-energy behavior.

    
% Q3 - Trivial Observer (Full State Measured)

    %   State Estimation Error:
    %       Compared to Q2, only the State Estimation Errors show important structural differences.
    %       The profiles show an immediate drop from the initial offset values directly to 0 at the very first time step, since we are no longer integrating a differential error equation.
    %       This behavior derives from the fact that the full state is measured; we are not predicting the missing variables anymore, but we are directly 'reading' them from the sensors (xHat = y).
    %       The final numerical error slightly increases due to the integration time step (dt=0.01) of the real state.
    
    %   Other Graphs:
    %       The State, Trajectory, and Input graphs do not change at all compared to Q2, because the optimal control input u(t) is generated a priori from matrices A, B and the final state x_f, 
    %       completely independent of the chosen observer architecture.