% Numerical Programming 2 CSE
% Tutorial 2 - ODE - Runge-Kutta and Multistep Methods
% Author: Bernhard Gatzhammer

% Matlab code for exercise 1c)

clear all
close all

y0 = [3/4*pi, 0]; % Initial angular displacement and angular velocity
tend = 10.0;      % End time of simulation
n = 1000;         % Number of timesteps
t = linspace(0, tend, n);
y = pendulum_rk(t, y0);
y = (y/(2*pi)) * 360; % Transforming angular radiant to angular degree
plot(t, y(:,1),'r') 
xlabel('t [sec]')
ylabel('displacement angle [?]')
grid on