% Numerical Programming 2 CSE
% Tutorial 2 - ODE - Runge-Kutta and Multistep Methods
% Author: Bernhard Gatzhammer

% Matlab code for exercise 1c)

function y = pendulum_rk(t, y0)
    %PENDULUM_RK Solves the system of ODEs describing the pendulum motion by 
    %the Runge-Kutta method
    %
    % Parameters:
    % t: Discrete times to compute solutions for.
    % y0: Vector of initial values [angle, angular velocity] at t(1)
    %
    % Return values:
    % y: Vector of integrated values (two row matrix)

    g = 9.81; % Gravitational constant
    l = 0.6;  % Length of pendulum
    f = @(dy) [-g/l*sin(dy(2)), dy(1)]; % Right-hand side of system of ODEs
    n = length(t);
    y = zeros(n, 2); % Allocate output function values
    y(1,:) = y0; 
    dt = t(2) - t(1);
    for k=2:n;
        % Compute stages
        f1 = f(y(k-1,:));
        f2 = f(y(k-1,:) + dt*0.5*f1);
        f3 = f(y(k-1,:) + dt*(-f1 + 2*f2));
        % Compute next timestep
        y(k,:) = y(k-1,:) + dt*(1/6 * f1 + 4/6 * f2 + 1/6 * f3);
    end
end

