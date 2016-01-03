% Numerical Programming 2 CSE
% Tutorial 12 - PDE - Finite Element Method 3
% Author: Bernhard Gatzhammer

% Matlab code to solve exercise 1b)

function ex1b()
    close all
    clear all

    i = 3;
    Nelem = 2^i;
    h = 1/Nelem;
    A = zeros(Nelem-1,Nelem-1);
    b = h * ones(Nelem-1,1);
    for i=1:Nelem-1
        A(i,i) = 2;
        if i > 1
            A(i-1,i) = -1;
        end
        if i < Nelem-1
            A(i+1,i) = -1;
        end
    end
    A = (1/h) * A;
    
    % Case with homogenous boundary conditions
    a1 = 0;
    b1 = 0;
    u = A\b;
    coords = linspace(0,1,Nelem+1);
    plot(coords', [a1; u; b1], 'b+-');
    hold on

    % Case with equal boundary conditions
    a2 = 1;
    b2 = 1;
    b = h * ones(Nelem-1,1);
    b(1) = b(1) + a2/h;
    b(end) = b(end) + b2/h;
    u0 = A\b;    
    % The boundary functions are added implicitely in the plot
    plot(coords', [a2; u0; b2], 'r+-');

    % Case with different boundary conditions
    a3 = 0;
    b3 = 1;
    b = h*ones(Nelem-1,1);
    b(1) = b(1) + a3/h;
    b(end) = b(end) + b3/h;
    u0 = A\b;
    % The boundary functions are added implicitely in the plot
    plot(coords', [a3; u0; b3], 'g+-');
    
    legend('homogenous', 'constant', 'linear')
    title('discrete boundary functions')
    grid on
end