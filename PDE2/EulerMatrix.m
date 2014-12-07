function matrix = EulerMatrix(Nx,Ny,dt)
% EULERMATRIX a function that returns
%   the matrix using the Euler Method
%   with dimensions Nx and Ny respectively.

b = (Nx+1)^2;
c = (Ny+1)^2;
alfaX= dt*b;  % alfaX= dT/ hx^2
alfaY=dt*c;    % dT/ hy^2
A = repmat(1-2*(alfaX+alfaY), 1, Nx);
B = repmat(alfaX, 1, Nx-1);
C = repmat(alfaY, 1, Nx*(Ny-1));

matrix = [];
temp = diag(A) + diag(B, -1) + diag(B, 1);
for i=1:Ny
    matrix = [matrix; zeros(Nx,(i-1)*Nx), temp, zeros(Nx, (Ny-i)*Nx)];
end
matrix = matrix + diag(C, -Nx) + diag(C, Nx);
    