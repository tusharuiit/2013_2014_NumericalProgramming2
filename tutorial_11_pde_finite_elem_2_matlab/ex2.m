% Numerical Programming 2 CSE
% Tutorial 11 - PDE - Finite Element Method 2
% Author: Bernhard Gatzhammer

% Matlab code to solve exercise 2

function ex2()
    i = 5;

    N = 2^i;
    h = 1/N;
    dofs = (N-1)*(N-1);
    A = zeros(dofs, dofs);
    b = zeros(dofs,1);
    % Reference element matrix part of left-hand side gradient integrals
    Aep = [ 2/3 -1/6 -1/6 -1/3;
           -1/6  2/3 -1/3 -1/6;
           -1/6 -1/3  2/3 -1/6;
           -1/3 -1/6 -1/6  2/3];
    % Reference element matrix part of left-hand side non-gradient integrals
    Aeq = [1/9  1/18 1/18 1/36;
           1/18 1/9  1/36 1/18;
           1/18 1/36 1/9  1/18;
           1/36 1/18 1/18 1/9 ];
    % Construction of global element matrix from reference element matrices
    Ae = Aep + h^2 * Aeq;
    % Construction of global stiffness matrix in loop over elements
    for j=1:N
        for i=1:N
            e = (i-1) + (j-1)*N;
            dof = e-N-j+2;
            eft = [dof, dof+1, dof+N-1, dof+N];
            
            % Mark boundary nodes
            if i==1
                eft(1) = -1;
                eft(3) = -1;
            end
            if i==N
                eft(2) = -1;
                eft(4) = -1;
            end
            if j==1
                eft(1) = -1;
                eft(2) = -1;
            end
            if j==N
                eft(3) = -1;
                eft(4) = -1;
            end
            
            % Compute rhs
            xi = (i-1)*h;
            yi = (j-1)*h;
            be = ones(1,4) * (xi^2 + yi^2)/4;
            be(1) = be(1) + h^2/12 + (yi/6 + xi/3)*h;
            be(2) = be(2) + h^2/6 + (yi/3 + xi/6)*h;
            be(3) = be(3) + h^2/6 + (yi + xi)*h;
            be(4) = be(4) + h^2/4 + (yi/3 + xi/3)*h;
            be = h^2 * be;
            
            % Add up element entries to global matrix and rhs
            for k=1:4
                for l=1:4                    
                    if eft(k) ~= -1 && eft(l) ~= -1
                        A(eft(k),eft(l)) = A(eft(k),eft(l)) + Ae(k,l);
                    end
                end
                if eft(k) ~= -1
                    b(eft(k)) = b(eft(k)) + be(k);
                end
            end
        end
    end
%     A
%     b
    u = A\b; % Solve system
    
    % plot results
    coords = linspace(0,1,N+1);
    [xmesh, ymesh] = meshgrid(coords);
    umesh = zeros(N+1,N+1);
    for i=1:N-1
        for j=1:N-1
            dof = i + (j-1)*(N-1);
            umesh(i+1,j+1) = u(dof);
        end
    end
    figure
    surf(xmesh, ymesh, umesh);
end