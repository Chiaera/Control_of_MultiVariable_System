%% Q4. Euler Estimation (step integration)

% Starting from the state equation x.(k) = Ax + Bu
% compute the next integration step: x(k+1) = x(k) + dt*x.(k)

function [x] = SimulateSysEuler(A, B, u, t, x0, dt)
    x = zeros(size(A,1),size(t,2)); % [6xN]
    x(:,1) = x0; % initial state

    % loop 
    for k = 1 : size(t,2) - 1
        x(:,k+1) = x(:,k) + dt*(A*x(:,k) + B*u(:,k));
    end
end


%% GENERAL OBSERVATIONS:
% - computational lighter than Analytical but accurate only if dt small enough
% - can diverge if the system is unstable or dt is too large