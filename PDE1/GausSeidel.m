function x = GausSeidel(b, Nx, Ny)
% GAUSSEIDEL performs Gaus-Seidel iteration
%   of the given system, returning matrix of
%   computed approximations at grid points.
% b .... right side of the system
% Nx ... x dimension of the system matrix
% Ny ... y dimension of the system matrix

% Initializing nonzero values of the matrix:
B = (Nx+1)^2;
C = (Ny+1)^2;
A = -2*(B+C);

% Setting starting approximation matrix:
x = zeros(Nx+2, Ny+2);

% Setting starting residual norm:
R = 1;
N = Nx*Ny;

while R>1e-4
    Rtemp = 0;
    for j = 2:Ny+2
        for i = 2:Nx+1
            if j<Ny+2
                x(i,j) = 1/A*(b((j-2)*Nx+i-1)-C*(x(i,j+1)+x(i,j-1))-B*(x(i+1,j)+x(i-1,j)));
            end
            if j>2
                Rtemp = Rtemp + (b((j-3)*Nx+i-1) - ...
                    (C*(x(i,j)+x(i,j-2))+B*(x(i+1,j-1)+x(i-1,j-1))+A*x(i,j-1)))^2;
            end
        end
    end
    R = sqrt(1/N * Rtemp);
end

end