function x = GausSeidel(Nx, Ny, old, dt)
% GAUSSEIDEL performs Gaus-Seidel iteration
%   of the given system (in our case for a
%   step of Euler) at certain point in 
%   time, returning matrix of computed 
%   approximations at grid points.
%
% Nx ... x dimension of the system matrix
% Ny ... y dimension of the system matrix
% old .. matrix of soution at previous time
% dt ... time step


% Initializing nonzero values of the matrix:
B = dt*(Nx+1)^2;
C = dt*(Ny+1)^2;
A = 1+2*dt*(B+C);

% Setting starting approximation matrix:
x = zeros(Nx+2, Ny+2);

% Setting starting residual norm:
R = 10;
N = Nx*Ny;

while R>1e-6
    Rtemp = 0;
    for j = 2:Ny+1
        for i = 2:Nx+1
            x(i,j) = 1/A*(old(i,j)+C*(x(i,j+1)+x(i,j-1))+B*(x(i+1,j)+x(i-1,j)));

        end
    end

    for j = 2:Ny+1
	for i = 2:Nx+1
            Rtemp = Rtemp + (old(i,j)+C*(x(i,j+1)+x(i,j-1))+B*(x(i+1,j)+x(i-1,j))-A*x(i,j))^2;
        end
    end

    R = sqrt(1/N * Rtemp);

end
