%   Worksheet 3
% (I. Tominec, J. C. Medina, E. Breznik)
%
% -Solving partial differential equations.
% -----------------------------------------------
clear all; close all;

% Set parameter b of the problem:
b = @(Nx, Ny, d) -2*pi^2*(sin(pi*repmat(1:Nx, 1, Ny)/(Nx+1)).*sin(pi*d(1:end)/(Ny+1)));

% b) Building the system matrix is done by SystemMatrix, e.g.
% Nx = 5;
% Ny = 5;
% M = SystemMatrix(Nx,Ny)

% c) Gaus-Seidel iterative solver call:
% grid = GausSeidel(b, Nx, Ny)

% d) Solving the given system in three ways. + e) Visualization of solutions.
N = [7, 15, 31, 63];

%We also calculate runtimes and storage for f):
runtimes = zeros(3,4);
storage = zeros(3,4);
%
%For the g. part of the worksheet, some error values are calculated here:
errors = zeros(2,5);
%
%plotting:
figure(1);
set(gcf,'numbertitle','off','name','Direct solver with full matrix');
figure(2);
set(gcf,'numbertitle','off','name','Direct solver with sparse matrix');
figure(3);
set(gcf,'numbertitle','off','name','Gauss-Seidel solver');

for i=1:4
    Nx = N(i);
    Ny = N(i);
    %this we need later for plotting:
    [x,y] = meshgrid(1/(Nx+1)*(0:(Nx+1)), 1/(Ny+1)*(0:(Ny+1)));
    %
    tempvec = repmat(1:Ny,Nx,1);  % d = repmat(1:Ny,Nx,1)
    be = b(Nx,Ny, tempvec)';

    M = SystemMatrix(Nx,Ny);
    m = sparse(M);
    
    %   1. by MATLAB direct solver:
    tic; sol1 = M\be; runtimes(1,i)=toc;
    storage(1,i) = (Nx*Ny)^2 + 2*Nx*Ny;
    %plot:
    sol1 = fliplr(reshape(sol1, [Ny,Nx]))';
    Temp = [zeros(1,Nx+2);zeros(Ny,1),sol1,zeros(Ny,1);zeros(1,Nx+2)];
    figure(1);
    subplot(2,4,i);
    surf(x,y,Temp);
    title({'Nx,Ny=',Nx});
    subplot(2,4,4+i);
    contour(Temp);
    
    %clear some memory: 
    clear sol1;
    clear Temp;
    
    %   2. by direct solver on a sparse matrix:
    tic; sol2 = m\be; runtimes(2,i)=toc;
    storage(2,i) = (5*Nx*Ny - 2*Ny - 2*Nx) + 2*Nx*Ny; % bcs we only store nonzero elements
    %plot:
    sol2 = fliplr(reshape(sol2, [Ny,Nx]))';
    Temp = [zeros(1,Nx+2);zeros(Ny,1),sol2,zeros(Ny,1);zeros(1,Nx+2)];
    figure(2);
    subplot(2,4,i);
    surf(x,y,Temp);
    title({'Nx,Ny=',Nx});
    subplot(2,4,4+i);
    contour(Temp);
    
    %clear some memory: 
    clear sol2;
    clear Temp;
    
    %   3. by Gaus-Seidel:
    tic; sol3 = GausSeidel(be, Nx, Ny); runtimes(3,i)=toc;
    storage(3,i)= (Nx+2)*(Ny+2) + Nx*Ny;
    %plot:
    figure(3);
    subplot(2,4,i);
    surf(x,y,sol3);
    title({'Nx,Ny=',Nx});
    subplot(2,4,4+i);
    contour(sol3);

    
    %%%%%%%%%% this part includes calculations for g)
    error = 0;
    for v=1:Ny
        for s=1:Nx
            error = error + (sol3(v+1,s+1) - sin(pi*v/(Ny+1))*sin(pi*s/(Nx+1)))^2;
        end
    end 
    errors(1,i) = sqrt(1/(Nx*Ny) * error);
    %%%%%%%%%%
    
    %clear some memory:
    clear sol3;

end


% f) Runtime and storage comparisons. TABLES:
f=figure;
tabgp = uitabgroup(f,'Position',[.05 .5 .95 .45]);
cnames = {'Nx,Ny=7','Nx,Ny=15','Nx,Ny=31','Nx,Ny=63'};
rnames = {'Runtime','Storage'};

tab1 = uitab(tabgp,'Title','FullMatrix');
data = [runtimes(1,1:end); storage(1, 1:end)];
% Create the uitable
t = uitable(tab1,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);


tab2 = uitab(tabgp,'Title','SparseMatrix');
data = [runtimes(2,1:end); storage(2, 1:end)];
% Create the uitable
t = uitable(tab2,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);


tab3 = uitab(tabgp,'Title','Gaus-Seidel');
data = [runtimes(3,1:end); storage(3, 1:end)];
% Create the uitable
t = uitable(tab3,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);



% g) Errors.
% (now it only remains to calculate for n=127, and error reductions)
Nx=127; Ny=127;
tempvec = repmat(1:Ny,Nx,1);
be = b(Nx,Ny, tempvec);

sol3 = GausSeidel(be, Nx, Ny);

error = 0;
for v=1:Ny
    for s=1:Nx
        error = error + (sol3(v+1,s+1) - sin(pi*v/(Ny+1))*sin(pi*s/(Nx+1)))^2;
    end
end 
errors(1,5) = sqrt(1/(Nx*Ny) * error);

%lets clean some memory:
clear sol3;

%error reduction factors:
%(when making table: errors(2,1)='-')
for k = 2:5
    errors(2,k) = errors(1,k-1)/errors(1,k);
end

%TABLES
f=figure;
tabgp = uitabgroup(f,'Position',[.05 .5 .95 .45]);
cnames = {'Nx,Ny=7','Nx,Ny=15','Nx,Ny=31','Nx,Ny=63', 'Nx,Ny=127'};
rnames = {'Error', 'Error red.'};

tab = uitab(tabgp,'Title','Gauss-Seidel Errors');
data = [errors(1,1:end); NaN, errors(2, 2:end)];
% Create the uitable
t = uitable(tab,'Data', data,'ColumnName', cnames,'RowName', rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);

    