%                WORKSHEET 3- Solving a PDE  (E. Breznik, J. C. Medina, I. Tominec)
%
% This program solves
close all; clear all;

% The problem definition.
T_der = @(x, y)(-2.*pi^2.*sin(pi.*x).*sin(pi.*y));

N=[7,15,31,63];
figure(1);
figure(2);
figure(3);
plot_counter=1;
for i=1:4
% b) Building the system matrix.
Nx = N(i);
Ny = N(i);
hx = 1/(Nx+1);
hy = 1/(Ny+1);
M = SystemMatrix(Nx,Ny);



% c) Gaus-Seidel iterative solver.


% d) Solving the given system in three ways.
%Finding the function vector of the second derivative

[ x y ] = meshgrid( hx*(0:(Nx+1)), hy*(0:(Ny+1)) );
f_der=T_der(x,y); %Evaluating the derivative
f_der=f_der(2:Nx+1,2:Ny+1);   %Matrix to column, taking out the zeros
f_der=f_der(:);

%-   1. by MATLAB direct solver:
tic

Temperature=M\f_der;
Temperature=vec2mat(Temperature,Nx); %Vector to column. Columns are determined only by the Nx
TempFull= [zeros(1,Nx+2);zeros(Ny,1),Temperature,zeros(Ny,1);zeros(1,Nx+2)]; %Inserting the boundary given conditions

time_Full(i)=toc;
%-   2. by direct solver on a sparse matrix:
tic

M=sparse(M);
Temperature=M\f_der;
Temperature=vec2mat(Temperature,Nx); %Vector to column. Columns are determined only by the Nx
TempSparse= [zeros(1,Nx+2);zeros(Ny,1),Temperature,zeros(Ny,1);zeros(1,Nx+2)]; %Inserting the boundary given conditions

time_Sparse(i)=toc;
%   3. by Gaus-Seidel:

%- e) Visualization of solutions.
figure(1);
subplot(2,4,plot_counter)
surf(x,y,TempFull)
subplot(2,4,plot_counter+1)
contour(TempFull)
title({'Nx,Ny=',Nx});

figure(2);
subplot(2,4,plot_counter)
surf(x,y,TempSparse)
subplot(2,4,plot_counter+1)
contour(TempSparse)
title({'Nx,Ny=',Nx});

plot_counter=plot_counter+2;

%- f) Runtime comparisons.



% g) Errors.
end

%TABLES

f=figure;
tabgp = uitabgroup(f,'Position',[.05 .5 .95 .45]);
cnames = {'Nx,Ny=7','Nx,Ny=15','Nx,Ny=31','Nx,Ny=63'};
rnames = {'Runtime','Storage'};

tab1 = uitab(tabgp,'Title','FullMatrix');
data = [time_Full; zeros(1,4)];
% Create the uitable
t = uitable(tab1,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);



tab2 = uitab(tabgp,'Title','SparseMatrix');
data = [time_Sparse; zeros(1,4)];
% Create the uitable
t = uitable(tab2,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);



