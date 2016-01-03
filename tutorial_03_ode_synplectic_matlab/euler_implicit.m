% Numerical Programming 2 CSE
% Tutorial 3 - ODE - Synplectic Methods
% Author: Bernhard Gatzhammer

% Matlab code for implicit Euler solving the hamiltonian ODE of exercise 1

function y = euler_implicit(odef, t, y0)
    y = zeros(length(t),2);
    y(1,:) = y0;
    dt = t(2) - t(1);
    for i=2:length(t)
        p_old = y(i-1,1);
        q_old = y(i-1,2);
        f = @(y) [y(1) - p_old - dt*cos(y(2)); y(2) - q_old - dt*y(1)];
        options = optimset('TolX', 1e-10, 'Display', 'off');
        y(i,:) = fsolve(f,y(i-1,:),options);
    end
end

