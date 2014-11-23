function poisson2d

% Define the mesh
N = 50;
M = 65;
a = 0; b = 2; c = -1 ; d = 2;
hx = (b-a) / (N+1); hy = (d-c) / (M+1);
[ y x ] = meshgrid( c + hy*(0:(M+1)), a + hx*(0:(N+1)) );

% Number of iterations
maxIter = 800;

%S equal to right hand sight of differential equation
% Define and plot source density S(x, y)
S = 5 * x.^2 .* y;
figure(2)
subplot(2, 2, 1);
imagesc(x(2:end-1, 1), y(1, 2:end-1), S')
title('S')

% Define the solution vector f as a matrix in Matlab.
% Insert Dirichlet BC into array, e.g. f = x on all boundaries
% We set it to x in the entire domain, interior points will be
% modified during the iteration.
f = x;


% Compute residual r
norm_r = residualNorm(f, S, hx, hy)

for loop = 1:maxIter
   % f = relaxJacobi(f, S, hx, hy);
    f = relaxGS(f, S, hx, hy);
    if(0 == mod(loop, 10))
        % Every 10 iterations, plot the solution and residual etc.
        subplot(2, 2, 2);
        imagesc(x(2:end-1, 1), y(1, 2:end-1), f')
        title('f')
        [ no, resi ] = residualNorm(f, S, hx, hy);
        norm_r(end+1) = no;
        no
        subplot(2, 2, 3);
        imagesc(x(2:end-1, 1), y(1, 2:end-1), resi')
        title('Residual');
        subplot(2, 2, 4);
        semilogy(norm_r);
        title('||Residual||');
        pause(0.2)
    end
end
end

function [ norm_r, r ] = residualNorm(f, S, hx, hy)
    % Vector of inner indices to access element of f, S
    J = 2:(size(f, 1)-1);
    K = 2:(size(f, 2)-1);
   
    r = S(J, K) - (f(J+1, K) - 2*f(J, K) + f(J-1, K) )/hx^2 ...
        - (f(J, K+1) - 2*f(J, K) + f(J, K-1) )/hy^2;
    norm_r = sqrt( sum(sum(r.^2)) / (length(J) * length(K)) );
end

function f = relaxJacobi(f, S, hx, hy)
    % return new value f of iteration l+1 from old value at l
    J = 2:(size(f, 1)-1);
    K = 2:(size(f, 2)-1);
   
      f(J, K) = (hy^2*(f(J-1,K) + f(J+1,K)) + hx^2*(f(J,K-1) + f(J,K+1))...
      -(hx*hy)^2*S(J,K))/(2*(hx^2+hy^2));
end

function f = relaxGS(f, S, hx, hy)
    J = 2:(size(f, 1)-1)
    K = 2:(size(f, 2)-1);
   
    for k = K
        for j = J
            f(j, k) = (hy^2*(f(j-1,k) + f(j+1,k)) + hx^2*(f(j,k-1) + f(j,k+1))...
      -(hx*hy)^2*S(j,k))/(2*(hx^2+hy^2));
 
        end
    end
end