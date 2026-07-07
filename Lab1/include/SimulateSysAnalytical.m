%% Q3. Analytical Estimation (differential equation)

% Calculates how the initial state evolves on its own (matrix exponential):
% - Adds, step by step, the effect of each past input 
% - Each input is "propagated" into the future thru the matrix exponential 
% - The sum of all these contributions gives the exact solution.

function [x] = SimulateSysAnalytical(A, B, u, t, x0, dt)
    % define trajectory matrix [6xN]
    x = zeros(size(A,1),size(t,2)); 
    x(:,1) = x0; 
    
    % for every k: x(k) = exp(At)x0+integ_k()
    for k = 2 : size(t,2)
        I_k = 0; % initialization integral term k [6x1]
        
        % Riemann sum to approximate the integral
        for i = 1 : k - 1
            I_k = I_k + expm(A*(t(k)-t(i)))*B*u(:,i)*dt;
        end
        x(:,k) = (expm(A*(t(k)))*x0 + I_k);
    end
end



%% GENERAL OBSERVATIONS:
% - directly uses the matrix exponential -> exact solution, 
% - calculates the effect of the input thru a convolution integral 