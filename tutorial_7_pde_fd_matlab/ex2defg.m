% Numerical Programming 2 CSE
% Tutorial 7 - PDE - Finite Difference Method
% Author: Bernhard Gatzhammer

% Matlab code for exercise 2d), e), f) and g)

close all
clear all

b = -1;
u0 = 0;
un = 1;
i_begin = 4;
i_end = 11;
errorsL = zeros(i_end-i_begin+1,1);
errorsLp = zeros(i_end-i_begin+1,1);
hs = zeros(i_end-i_begin+1,1);
for i=i_begin:i_end
    h = 1/2^i;
    n = floor(1/h + 0.5) - 1;
    L = zeros(n,n);
    Lp = zeros(n,n);
    u = zeros(n,1);
    
    % Create matrix entries
    Lleft = -(1/h^2 + b/(2*h));
    Lmid = 2/h^2;
    Lright = -(1/h^2 - b/(2*h));
    Lpleft = -1/h^2;
    Lpmid = 2/h^2 - b/h;
    Lpright = -(1/h^2 - b/h);
    
    % Set matrix entries
    L(1,1) = Lmid;
    L(1,2) = Lright;
    Lp(1,1) = Lpmid;
    Lp(1,2) = Lpright;
    for l=2:n-1
       L(l,l-1) = Lleft;
       L(l,l)   = Lmid;
       L(l,l+1) = Lright;
       Lp(l,l-1) = Lpleft;
       Lp(l,l)   = Lpmid;
       Lp(l,l+1) = Lpright;
    end
    L(n,n-1) = Lleft;
    L(n,n) = Lmid;
    Lp(n,n-1) = Lpleft;
    Lp(n,n) = Lpmid;
    
    % Create and set right-hand sides
    bL = zeros(n,1);
    bLp = zeros(n,1);
    bL(1) = (1/h^2 + b/(2*h))*u0;
    bL(n) = (1/h^2 - b/(2*h))*un;
    bLp(1) = (1/h^2)*u0;
    bLp(n) = (1/h^2 - b/h)*un;
    
    % Solve systems
    uL = L\bL;
    uLp = Lp\bLp;
    
    % Compute analytical solution at discrete points
    factor = 1/(1-exp(b));
    for l=1:n
        u(l) = factor * (1 - exp(b*l*h));
    end
    
    % Compute errors
    errorsL(i-i_begin+1) = max(abs(u - uL));
    errorsLp(i-i_begin+1) = max(abs(u - uLp));    
    hs(i-i_begin+1) = h;
    
    % Exercise 2g), set b=-1e4 also
    if i == 5
        ns = transpose(linspace(0+h,1-h,n));
        plot(ns, uL, 'b-')
        hold on
        plot(ns, uLp, 'r-')
        plot(ns, u, 'g-')
        xlabel('h')
        ylabel('u')
        legend('uL', 'uL+', 'u')
    end
end

figure
loglog(hs, errorsL, 'b+-')
hold on
loglog(hs, errorsLp, 'r+-')
xlabel('h')
ylabel('max error')
legend('L', 'L+')
grid on