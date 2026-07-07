%% MVCS LAB2
clc
close all
clear

%% Set SCENARIO 
% 1 = Q1.5 - 4 thrusters, redundant controllable system
% 2 = Q1.6 - No thruster 2, system still controllable
% 3 = Q1.y - No thruster 2 and 3, UNcontrollable system

scenario = 3;

%% Define matrices
% Q1.1 - Define the system
alfa = 1.4;
beta = 1.8;
gamma = 1.05;
m = 13.5;
J = 0.37;

% considering the state x=[x y theta vx vy omega]^T 

A = [ 0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      0 0 0 -alfa/m 0 0;
      0 0 0 0 -beta/m 0;
      0 0 0 0 0 -gamma/J ];

TAM = [ 0.7071 0.7071 -0.7071 -0.7071;
        -0.7071 0.7071 -0.7071 0.7071;
        -0.1888 0.1888 0.1888 -0.1888 ];

% Correct scenario configuration
switch scenario
    case 1
        disp('Q1.5');
    case 2 % Q1.6 - thruster 2 disables
        disp('Q1.6');
        TAM(:,2) = 0; 
    case 3 % Q1.7 - thruster 3 disables too
        disp('Q1.7');
        TAM(:,2) = 0;
        TAM(:,3) = 0; 
end

B_force = [ 0 0 0;   % input matrix wrt the force
            0 0 0; 
            0 0 0; 
            1/m 0 0; 
            0 1/m 0;
            0 0 1/J ];

B = 30 * B_force * TAM; % input matrix considering the TAM [6x4]


%%  Q1.2 - Controllability

% Compute Gramian matrix G(0,t) = integral 0->t (e^(A(t-tao) * BB^T * e^(A^T(t-tao))
T = 20; % simulation time
dt = 0.01; % simulation step
xf = [-1 2 3 0 0 0]';
syms tau real
eA = expm(A*(T-tau)); % exponential matrix
dG = eA * B * B' * (eA)'; % integral part
G = eval(int(dG,tau,0,T));

% check the rank
rankG = rank(G);
disp(['rank(G) = ' num2str(rankG)]);

if rankG == size(A,1) % IF rank(G) == 6 -> full rank -> invertible Gramian
    iG = inv(G);  % Gramian inverse
    disp("G is full rank! System is completely CONTROLLABLE")
else    % IF rank(g) < 6 -> Gramian not invertible
    iG = pinv(G);     % Moore-Penrose pseudo-inverse
    disp("System is not controllable")

    % Q1.8 - SVD analysis for non controllable system 
    [U,D,V] = svd(G); 

    disp("Matrix D (singular values on diagonal):")
    disp(diag(D)); % Read diagonal as a column

    disp("Matrix U (space state directions):")
    disp(U);

    disp("Matrix V (space state input(motors)):")
    disp(V);
end

    
%% Q1.3 - Compute output to bring x_0 in x_f
    
% Define time vector
time = [0:dt:T]
time = time';
LL = length(time); % samples number

% Compute system input u = B^T * e^(A^T(T-t) * G(0,T)^1 * x_f)
u(tau) = B'*eA'*iG*xf;
u_time = u(time);
uvalue = [eval(u_time{1}), eval(u_time{2}), eval(u_time{3}), eval(u_time{4})]; 

%% Q4. Simulate system

% Initialize state
x0 = [0 0 0 0 0 0]'; % we start from the origin
x = zeros(6,LL); % this matrix will store the state evolution
x(:,1) = x0; % firt instant

% Simulate system (you can use Euler integration)
for k = 1:LL-1
    u_k = uvalue(k,:)'; % extract the col of the tensions at t=k
    xdot = A * x(:,k) + B * u_k; % derivate at state k (instant variations)
    x(:,k+1) = x(:,k) + dt * xdot; % Euler integration to find x(k+1) 
end

%% Q1.5 - Plot state evolution
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Set figure to full screen

% Plot inputs
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

disp("final state - ref state = " + norm(x(:,end) - xf)); % should be small
% results:
    % Q1.5 = 0.0024227
    % Q1.6 = 0.002623
    % q1.7 = 0.11967   


%% PLOT COMMENTS:

% Q1.5 - redundant motor
    % The system is full-rank, so it is fully controllable.

    % POSITION
    % x: Starts from 0 and touches the reference exactly at 20s. 
    %    The start is slow, which is possible due to inertia and friction.
    % y: Same reasons as x.
    % yaw: The moment of inertia is low (J=0.37) with respect to the mass (m=13.5), 
    %      so the velocity of rotation can be considered constant; because of this, the line is linear.
    
    % VELOCITY
    % x: Correct behavior of an input to achieve the fastest path in a given interval.
    % y: Same as x, but the maximal value is almost double compared to x. This happens because 
    %    the goal distance for y (2m) is bigger than the x one (-1m).
    % yaw rate: Considering a lower inertia, a little impulse is enough to move the ROV, 
    %           so we can have a smoother curve.
     

% Q1.6
    % The system remains full-rank, so it is fully controllable.

    % The thrusters are 4, but in this case, we control only 3 DoF. 
    % However, we can cover all the 6 states anyways. This means 3 motors are sufficient to have a controllable system.

    % What changes with respect to Q1.5 is the symmetry of the force; with only 3 motors, 
    % the force necessary to rotate and translate is different between the right and left parts of the ROV.
    % This is easily visible in the yaw rate state.


% Q1.7 - Q1.8
    % In this case, the system has rank(G) = 5. The rank is not full, so the system is not controllable.

    % From the SVD, it is possible to see that there is one uncontrollable direction.
    % Indeed, the last value of D is zero.
    % Considering the U and V matrices, column 6 shows the same values for vx and vy, meaning these velocities are stuck together.

    % Proved graphically, the plots do not respect the proper behavior seen in Q1.5 and Q1.6.
    % The Moore-Penrose pseudo-inverse permits the ROV to touch the reference at 20s, but the behavior is still uncontrollable.
