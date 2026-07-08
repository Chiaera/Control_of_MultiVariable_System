clear;
clc;
close all;

%% SET SIMULATIOM TYPE

% Scenario
% Q1, Q2: enableClosedLoop = true, enableClosedLoopWithObserver = false, fullStateMeasured = false
% Q3: enableClosedLoop = true, enableClosedLoopWithObserver = true, fullStateMeasured = false
% Q4: enableClosedLoop = true, enableClosedLoopWithObserver = true, fullStateMeasured = true

enableClosedLoop = true; 
enableClosedLoopWithObserver = true;
fullStateMeasured = true;

if enableClosedLoopWithObserver && fullStateMeasured
    desc = "closed loop with observer (Full State Measured - Q4)";
elseif enableClosedLoopWithObserver
    desc = "closed loop with observer (Q3)";
elseif enableClosedLoop
    desc = "closed loop (Q1, Q2)";
else
    desc = "open loop";
end

%% SET SYSTEM PARAMETERS
alfa = 1.4;
beta = 1.8;
gamma = 1.05;
m = 13.5;
J = 0.37;
A = [0,0,0,1,0,0;
    0,0,0,0,1,0;
    0,0,0,0,0,1;
    0,0,0,-alfa/m,0,0;
    0,0,0,0,-beta/m,0;
    0,0,0,0,0,-gamma/J];

xf = [-1;2;3;0;0;0];
testdesc = "Test";

B = [ 0 0 0; 
      0 0 0;
      0 0 0;
      1/m 0 0;
      0 1/m 0;
      0 0 1/J ];
usize = min(size(B)); % == 3 input

T = 20;
dt = 0.01;

% Output measurement matrix C based on Exercise
if fullStateMeasured
    C = eye(6); % Q4. All state misured
else
    C = zeros(3,6); % Q2, Q3. Only pose variables are measured 
    C(1:3,1:3) = diag([1,1,1]);
end

% Optimal OpenLoop input 
syms tau real

tic
eA = expm(A*(T-tau)); % transition state matrix
dG = eA*B*B'*eA'; % integral part of the Gramian
G = vpa(int(dG,tau,0,T)); % complete Gramian

rankG = rank(G);
disp(['rank(G) = ' num2str(rankG)]);

if rankG == 6  % full-rank
    disp("System has full-rank: fully CONTROLLABLE")
    iG = inv(G);
else
    iG = pinv(G);
    warning("System has rank " + rankG + "! It is not fully controllable.")
end

u(tau) = B'*eA'*iG*xf; % nominal optimal input (lower control energy)

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

%% SET CONTROLLER PARAMETERS
% Q1. Closed Loop assuming that state is fully measurable

% Set Kc (controlability) - to stabilize ROV around target
Kc = zeros(3,6);

% Set parameters
Tau = 1;
xi = 1;
wn = 1;

% controllor gain
kc1 = -m * wn^2;
kc2 = -m * wn^2;
kc3 = -J * wn^2;
kc4 = alfa - 2*xi*wn*m;
kc5 = beta - 2*xi*wn*m;
kc6 = gamma - 2*xi*wn*J;

% Kc matrix
Kc(1:3,1:3) = diag([kc1, kc2, kc3]);
Kc(1:3,4:6) = diag([kc4, kc5, kc6]);

%% SET OBSERVER PARAMETERS
% Q2. Luenberger Observer assuming that pose variables only are measured

% Set Ko (observability) - to make the estimation converge

Tau = 1;
xi = 1;
wn = 1;
Ko = zeros(6,3);
wn = 10; % compared to control matric Kc: increased speed

% observer gain
kx1 = 2*xi*wn - alfa/m;
kx2 = 2*xi*wn - beta/m;
kx3 = 2*xi*wn - gamma/J;
kx4 = wn^2 - kx1*alfa/m;
kx5 = wn^2 - kx2*beta/m;
kx6 = wn^2 - kx3*gamma/J;

% Ko matrix based on Exercise
if fullStateMeasured
    % Q4. Closed Loop and Luenberger assuming that the whole state is measured 
    % xHat = y => Ko is multiplier-free
    Ko_loop = zeros(6,6); 
else
    Ko_loop = zeros(6,3);
    Ko_loop(1:3,1:3) = diag([kx1, kx2, kx3]);
    Ko_loop(4:6,1:3) = diag([kx4, kx5, kx6]);
end

%% SIMULATION
u_stdev = 0.1;

% sensor noise deviation based on Exercise
if fullStateMeasured
    y_stdev = 0.05;
else
    y_stdev = 0.0;
end

x_ol = x * 0; % ol means Open Loop
x_ol(:,1) = x0;

xHat(:,1) = x0 + [1, -1, 1.5, 0.5, -0.5, -1.5]'; % Initial estimate error

y = zeros(min(size(C)),LL);
yHat = zeros(min(size(C)),LL);
y_noise = zeros(min(size(C)),LL);
yHat(:,1) = C * xHat(:,1);

for k = 1:LL-1
    %% Control

    % OpenLoop: ROV follows optimal u(t) pre-computed
    % ClosedLoop: feedbak of real state x to correct the trajectory
    % ClosedLoop + Observer: compute feedback with estimate xHat

    u_k = uvalue(k,:)'; 
    xdot_ol_k = A*x_ol(:,k) + B * u_k;
    x_ol(:,k+1) = x_ol(:,k) + xdot_ol_k * dt; % robot state in OpenLoop
    
    if ~enableClosedLoop % Q1, Q2
        x(:,k+1) = x_ol(:,k+1);
        uref(:,k) = u_k; % OpenLoop considers only nominal u
    else
        if enableClosedLoopWithObserver % Q3 considering state estimation xHat
            uref(:,k) = u_k - Kc * (x_ol(:,k) - xHat(:,k)); % corrective feedback: choose between the real one x and the estimate one xHat
        else % Q1 real feedback of x
            uref(:,k) = u_k - Kc * (x_ol(:,k) - x(:,k));
        end 
        
        xdot_k = A * x(:,k) + B * uref(:,k); % real evolution of ROV in ClosedLoop
        x(:,k+1) = x(:,k) + xdot_k * dt; % Euler integration
    end

    %% Measurement update
    y(:,k) = C * x(:,k); 

    %% Observer
    if fullStateMeasured  % Q4 noisy readings for all 6 states 
        y_noise(:,k) = y(:,k) + y_stdev*randn(6,1);
       
        % Trivial estimation update for full measurement
        xHat(:,k+1) = y_noise(:,k);  
        yHat(:,k+1) = y_noise(:,k);
    else  % Q2,Q3 noisy readings for positions only 
        y_noise(1:2,k) = y(1:2,k) + y_stdev*randn(2,1);
        y_noise(3,k) = y(3,k) + y_stdev*randn;
        u_obs_k = uref(:,k);
    
        % Luenberg observer
        xdotHat = A * xHat(:,k) + B * u_obs_k + Ko_loop * (y_noise(:,k) - yHat(:,k));
        xHat(:,k+1) = xHat(:,k) + xdotHat * dt;
    
        % updated estimate output
        yHat(:,k+1) = C * xHat(:,k+1);
    end
end

%% PLOTS
% Plot states
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Set figure to full screen
sgtitle("State" + newline + desc) 
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
sgtitle("State estimation error" + newline + desc) 
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
legend('u_1','u_2','u_3','u_4')
xlabel('Time [s]')
ylabel('Inputs')
title(['Input ',newline,desc])
grid on;

% Assuming you have time vector t corresponding to your data
figure;
scatter(x(1,:), x(2,:), [], time, 'filled'); % Plot with color mapping to time
colormap(jet); % Set colormap
cb = colorbar; % Add colorbar
cb.Label.String = 'time [s]'; % Add label to colorbar
xlabel('x [m]');
ylabel('y [m]');
title(['Trajectory', newline, desc]); % Corrected title
grid on;

disp("Final state error:")
x(:,end) - xf

disp("Final est error:")
xHat(end,k+1) - x(end,k+1)


%% PLOTS COMMENTS

% GENERAL - comparison Lab 3 with Lab 4
% The nominal open-loop optimal input computed via the Gramian is exactly the same for all cases (Q1, Q3, Q4) because the system matrices and T are identical. 
% The real difference in Lab 4 is the closed-loop feedback correction: -Kc*(error). 
% The real trajectory changes based on the information used to calculate this error (true state, estimated state, or noisy state).


% Q1 - Q2: Closed loop (Real State) & Luenberger Observer
    % State & Trajectory: 
    %   The position profiles start from zero and smoothly touch the reference at T=20s.
    %   The velocity profiles touch the maximum at the middle and return to zero.
    %   The trajectory is a clear straight line, showing the minimum energetical path.
    %   These correct behaviours show the efficacy of Kc that corrects the nominal optimal trajectory at every instant using perfect true state data.
    
    % Input: 
    %       All the tensions are smooth and continue, without disconnection, confirming the optimal energy used.
    
    % State Estimation Error (Q2):
    %   We assume to measure only pose variables. The observer starts with a high given initial error. 
    %   Because we set the observer much faster than the controller (wn=10 vs wn=1), it converges to zero very fast (around 0.4s). 
    %   The initial spikes in the velocity errors represent the "peaking phenomenon", caused by the high gain of the observer reacting strongly to the initial offset.
    %   
    %   Final state error is on the order of 10^-3 so it is low. 
    %   Final est error is 0.


% Q3: Closed loop with Luenberger Observer (Estimated State)
    %   We don't measure velocities, so the feedback uses the estimated state xHat.

    % State and Trajectory (visible changes compared to Q1): 
    %   At the beginning, the observer has a huge initial error. Because the controller trusts this wrong xHat (give error), 
    %   the real state and trajectory show a small initial deviation (transient). 
    %   Since the observer converges very fast, the controller stabilizes quickly and reaches the target perfectly.

    %   This proves that we can design the controller and observer independently. 
    %   The closed loop works and is stable even starting from a bad initial estimation.

    %   Final state error is still on the order of 10^-3, just slightly different from Q1 due to the different initial dynamic transient.


% Q4 - Closed loop with Observer (Full State Measured + Noise)
    %   We assume the whole state is measured (C = eye(6)) but sensors have noise (y_stdev = 0.05). 
    %   The observer is trivial (xHat = y_noise).

    %   State and Input: 
    %   The position profiles successfully reach the target because the integration acts as a natural low-pass filter. 
    %   he velocities and the inputs are very noisy. 
    %   This happens because the controller reacts continuously to the raw noisy sensor data.

    %   Trajectory: 
    %       The path is mostly straight but shows micro-ripples caused by the noisy inputs continuously exciting the thrusters.

    %   State Estimation Error: 
    %       Unlike Q1 and Q3, the estimation error plot does not converge to a flat zero line, but forms a noise band centered at zero, 
    %       showing the stochastic noise injected at each step. Consequently, the printed "Final est error" is not exactly 0, 
    %       but a random value on the order of 10^-2 (matching the noise standard deviation), which perfectly simulates real-world noisy sensor readings.    

    