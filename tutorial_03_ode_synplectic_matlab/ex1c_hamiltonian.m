% Numerical Programming 2 CSE
% Tutorial 3 - ODE - Synplectic Methods
% Author: Bernhard Gatzhammer

% Matlab code for exercise 1c)

clear all
%close all

odef1 = @(t,y) cos(y(2)); % First ODE function handle
odef2 = @(t,y) y(1);      % Second ODE function handle
odef = @(t,y) [odef1(t,y); odef2(t,y)]; % System of first and second ODE f.h.
npoints = 100;
starttime = 0.2*pi;
endtime = starttime + (1*(8-starttime));
points0 = zeros(npoints,2);
points = zeros(npoints,2);
for n=1:npoints
    angle = (2*pi / npoints)*n;
    y0 = 0.4*[sin(angle), cos(angle)]; % Initial values
    points0(n,:) = y0; 
    t = linspace(starttime, endtime, 100);

    [t, y] = ode45(odef, t, y0); % Matlab explicit ODE solver
    %y = euler_explicit(odef, t, y0);
    %y = euler_synplectic(odef1, odef2, t, y0);
    
    points(n,:) = y(end,:);
end

plot(points0(:,1), points0(:,2),'b')
hold on
plot(points(:,1), points(:,2), 'r')
xlabel('p')
ylabel('q')
grid on