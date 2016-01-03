% Numerical Programming 2 CSE
% Tutorial 3 - ODE - Synplectic Methods
% Author: Bernhard Gatzhammer

% Matlab code for synplectic Euler variant 1 (from lecture) solving
% the hamiltonian ODE of exercise 1

function y = euler_synplectic(odef1, odef2, t, y0)
    y = zeros(length(t),2);
    y(1,:) = y0;
    dt = t(2) - t(1);
    for i=2:length(t)    
        y(i,2) = y(i-1,2) + dt*odef2(t(i), y(i-1,:));
        y(i,1) = y(i-1,1) + dt*odef1(t(i), [y(i-1,1), y(i,2)]);
    end
end