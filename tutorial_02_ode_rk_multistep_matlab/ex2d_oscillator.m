% Numerical Programming 2 CSE
% Tutorial 2 - ODE - Runge-Kutta and Multistep Methods
% Author: Bernhard Gatzhammer

% Matlab code for exercise 2d)

clear all
close all

y0 = [1, 0]; % Initial displacement and velocity

t = linspace(0, 2*pi, floor(2*pi*100));
y = oscillator_bdf2(t, y0);
y_exact = cos(t);
plot(t, y(:,1),'b')
grid on
xlabel('t [sec]')
ylabel('x(t) [m]')
figure
plot(t, y(:,1)'-y_exact, 'r')
xlabel('t [sec]')
ylabel('error(t) [m]')
grid on

% Perform convergence series and plot error
dts = [];
errors = [];
t_measure = 7/4 * pi;
for k=1:8
    dt = 1/(2^k);
    t = linspace(0, 5, 5*2^k);
    y = oscillator_bdf2(t, y0);
    dts = [dts; dt];
    n_measure = length(t);
    errors = [errors; abs(cos(t(n_measure)) - y(n_measure,1))];
end
figure
loglog(dts, errors, 'b+-')
xlabel('dt')
ylabel('error')
grid on
