% Numerical Programming 2 CSE
% Tutorial 3 - ODE - Synplectic Methods
% Author: Bernhard Gatzhammer

% Matlab code for exercise 1b)

clear all
close all

y0 = [0, 1];              % Initial values [p0, q0]
t = linspace(0, 20, 100); % Discrete times used in integration
Ht = zeros(length(t),1);  % Energy over time
odef1 = @(t,y) cos(y(2)); % First ODE function handle
odef2 = @(t,y) y(1);      % Second ODE function handle
odef = @(t,y) [odef1(t,y); odef2(t,y)]; % System of first and second ODE f.h.
H = @(t,y) 0.5*y(1)^2 - sin(y(2)); % Hamiltonian function handle
Ht(1) = H(t(1),y0);       % Initial energy measure
Ht(1)


[t, y] = ode45(odef, t, y0); % Matlab explicit ODE solver
%y = euler_explicit(odef, t, y0);
%y = euler_implicit(odef, t, y0);
%y = euler_synplectic(odef1, odef2, t, y0);

% Plot the oscillation
plot(t, y(:,1),'b') 
xlabel('t')
ylabel('x')
grid on

% Plot the energy error compared to the initial energy
for i=2:length(t)
    Ht(i) = H(t(i),y(i,:));
end
figure
plot(t,Ht - Ht(1),'r')
grid on
xlabel('t')
ylabel('energy error')
