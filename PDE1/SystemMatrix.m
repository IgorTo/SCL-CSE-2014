function matrix = SystemMatrix(Nx,Ny)
% SYSTEMMATRIX a function that returns
%   the matrix of the discussed system
%   with dimensions Nx and Ny respectively.

b = (Nx+1)^2;
c = (Ny+1)^2;
A = repmat(-2*(c+b), 1, Nx);
B = repmat(b, 1, Nx-1);
C = repmat(c, 1, Nx*(Ny-1));

matrix = [];
temp = diag(A) + diag(B, -1) + diag(B, 1);
for i=1:Ny
    matrix = [matrix; zeros(Nx,(i-1)*Nx), temp, zeros(Nx, (Ny-i)*Nx)];
end
matrix = matrix + diag(C, -Nx) + diag(C, Nx);
    