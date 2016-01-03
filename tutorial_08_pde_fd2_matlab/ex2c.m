% Numerical Programming 2 CSE
% Tutorial 8 - PDE - Finite Difference Method 2
% Author: Bernhard Gatzhammer

% Matlab code for exercise 2c)

close all 
clear all

kend = 6;
maxerrors = zeros(kend,1);
hs = zeros(kend,1);

sol = @(x,y) exp(pi*x)*sin(pi*y)+0.5*(x*y)^2;

for k=2:kend
    n = 2^k;
    h = 1/(n+1);
    L = zeros(n^2 ,n^2);
    x = zeros(n^2,1);
    exact = zeros(n^2,1);
    b = zeros(n^2,1);
   
    % Five-point star   
    for i=1:n^2    
        x = (mod(i-1,n) +1) * h;
        y = (floor((i-1)/n) +1) * h;
                
        % Initialize right-hand side
        b(i) = 3*(h^2)*(x^2 + y^2);
        
        % Lower boundary     
        if i < n+1
            b(i) = b(i) - 4*sol(x,y-h) + (1/4)*sol(x,y-2*h);
        elseif i < (2*n)+1
            b(i) = b(i) + 1/4*sol(x,y-2*h);
            L(i,i-n) = 4;
        else
            L(i,i-n) = 4;
            L(i,i-2*n) = -1/4;
        end
        % Upper boundary
        if i > n*(n-1)     
            b(i) = b(i) - 4*sol(x,y+h) + (1/4)*sol(x,y+2*h);
        elseif i > n*(n-2)      
            b(i) = b(i) + 1/4*sol(x,y+2*h);
            L(i,i+n) = 4;
        else
            L(i,i+2*n) = -1/4;
            L(i,i+n) = 4;
        end
        % Left boundary
        if mod(i-1,n)+1 < 2
            b(i) = b(i) - 4*sol(x-h,y) + (1/4)*sol(x-2*h,y);
        elseif mod(i-1,n)+1 < 3
            b(i) = b(i) + (1/4)*sol(x-2*h,y);
            L(i,i-1) = 4;
        else
            L(i,i-1) = 4;
            L(i,i-2) = -1/4;
        end
        % Right boundary        
        if mod(i-1,n) + 1 > n-1
            b(i) = b(i) - 4*sol(x+h,y) + (1/4)*sol(x+2*h,y);
        elseif mod(i-1,n) + 1 > n-2
            b(i) = b(i) + (1/4)*sol(x+2*h,y);
            L(i,i+1) = 4;
        else
            L(i,i+1) = 4;
            L(i,i+2) = -1/4;
        end
        
        % main diagonal
        L(i,i) = -15;        
        
        exact(i) = sol(x,y);
    end
    
    % Solve system
    x = L\b;
    
    % Store h and maximal error
    hs(k) = h;
    maxerrors(k) = max(abs(exact - x));
    
    if k == 2  % Perform solution plots
        coords = linspace(0,1,n+2);    
        [xx, yy] = meshgrid(coords, coords);
        xgrid = zeros(n+2,n+2);
        exactgrid = zeros(n+2,n+2);
        errors = zeros(n+2,n+2);
        for j=1:n+2
            for i=1:n+2
                exactgrid(j,i) = exp(pi*(i-1)*h)*sin(pi*(j-1)*h)+0.5*((i-1)*(j-1)*h^2)^2;
                xgrid(j,i) = exactgrid(j,i);
            end
        end
        for j=2:n+1
            for i=2:n+1
                index = i - 1 + (j-2)*n;
                xgrid(j,i) = x(index);
                %exactgrid(i,j) = exp(pi*(i-1)*h)*sin(pi*(j-1)*h)+0.5*((i-1)*(j-1)*h^2)^2;
                errors(j,i) = exactgrid(j,i) - x(index);
            end
        end
        surf(xx,yy,exactgrid);
        title('exact solution')
        xlabel('x')
        ylabel('y')
        figure
        surf(xx,yy,xgrid);
        title('numerical solution')
        figure
        surf(xx,yy,errors);
        title('error')
    end
end


% Plot maximal errors in dependence of h
figure
loglog(hs, maxerrors, 'b+-')
xlabel('h')
ylabel('max |error|')
grid on