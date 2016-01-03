% Numerical Programming 2 CSE
% Tutorial 8 - PDE - Finite Difference Method 2
% Author: Bernhard Gatzhammer

% Matlab code for exercise 1c)

close all
clear all

f = @(x,y) x*y*exp(x*y);
fxx = @(x,y) 2*y^2*exp(x*y) + x*y^3*exp(x*y);
fyy = @(x,y) 2*x^2*exp(x*y) + x^3*y*exp(x*y);
s1 = @(x,y,h) f(x-h,y) + f(x+h,y) + f(x,y-h) + f(x,y+h);
s2 = @(x,y,h) f(x-2*h,y) + f(x+2*h,y) + f(x,y-2*h) + f(x,y+2*h);
fd  = @(x,y,h) (4*s1(x,y,h) - (1/4)*s2(x,y,h) - 15*f(x,y))/(3*h^2);
x = 1;
y = 1;
exact = fxx(x,y) + fyy(x,y)
kend = 20;
numeric = zeros(1,kend);
hrange  = zeros(1,kend);
for k=1:kend
    hrange(k) = 1/2^k;
    numeric(k) = fd(x,y,hrange(k));
end
plot(hrange, numeric, '+-b')
xlabel('h')
ylabel('finite difference')
figure
error = abs(numeric - exact);
%plot(hrange, error, '+-r')
%semilogy(hrange, error, '+-r')
loglog(hrange, error, 'r+-')
%grid on
xlabel('h')
ylabel('|error|')