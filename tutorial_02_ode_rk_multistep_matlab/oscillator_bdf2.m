% Numerical Programming 2 CSE
% Tutorial 2 - ODE - Runge-Kutta and Multistep Methods
% Author: Bernhard Gatzhammer

% Solution to exercise 2d)

function y = oscillator_bdf2(t, y0)
    %OSCILLATOR_BDF2 Solves the system of ODEs describing the harmonic oscillator 
    %by the two-step BDF method
    %
    % Parameters:
    % t: Discrete times to compute solutions for.
    % y0: Vector of initial values [discplament, velocity]
    %
    % Return values:
    % y: Vector of integrated values (two column matrix)

    y = zeros(length(t), 2);

    % Set proper initial conditions for both old steps from exact solution
    dt = t(2) - t(1);
    y(1,:) = y0;
    y(2,:) = [y0(1)*cos(dt), -y0(2)*sin(dt)];

    % Time integrator system matrix
    A = [1      -dt*2/3;
         dt*2/3       1];

    for k=3:length(t);
        % Compute system rhs
        rhs = 4/3*y(k-1,:) - 1/3*y(k-2,:);
        % Solve system
        y(k,:) = (A\(rhs'))';
    end
end

