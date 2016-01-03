% Numerical Programming 2 CSE
% Tutorial 7 - PDE - Finite Difference Method
% Author: Bernhard Gatzhammer

% Matlab code for exercise 2b)

close all
clear all

b = -1;
i_begin = 4;
i_end = 11;
inftynormsL = zeros(i_end-i_begin+1,1);
inftynormsLp = zeros(i_end-i_begin+1,1);
hs = zeros(i_end-i_begin+1,1);
for i=i_begin:i_end
    h = 1/2^i;
    n = floor(1/h + 0.5) - 1;
    L = zeros(n,n);
    Lp = zeros(n,n);
    
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
    
    % Compute inverse and infinity norm
    Linv = L\eye(n,n);
    Lpinv = Lp\eye(n,n);
    rowsumsL = zeros(n, 1);
    rowsumsLp = zeros(n, 1);
    for l=1:n
        for k=1:n
            rowsumsL(l) = rowsumsL(l) + abs(Linv(l,k));
            rowsumsLp(l) = rowsumsLp(l) + abs(Lpinv(l,k));
        end
    end
    inftynormsL(i-i_begin+1) = max(rowsumsL);
    inftynormsLp(i-i_begin+1) = max(rowsumsLp);
    hs(i-i_begin+1) = h;
end

plot(hs, inftynormsL, 'b+-')
hold on
plot(hs, inftynormsLp, 'r+-')
xlabel('h')
ylabel('inftynorm')
legend('L', 'L+')