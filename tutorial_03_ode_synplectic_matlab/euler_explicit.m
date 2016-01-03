% Numerical Programming 2 CSE
% Tutorial 3 - ODE - Synplectic Methods
% Author: Bernhard Gatzhammer

% Matlab code for explicit Euler solving the hamiltonian ODE of exercise 1

function y = euler_explicit(odef, t, y0)
    y = zeros(length(t),2);
    y(1,:) = y0;
    dt = t(2) - t(1);
    for i=2:length(t)
        y(i,:) = y(i-1,:) + dt*transpose(odef(t(i), y(i-1,:)));
    end
end

