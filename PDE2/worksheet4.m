%Tables 
%Stability
stable=[1 1 1 1 1 1 1;
        0 0 1 1 1 1 1;
        0 0 0 0 1 1 1;
        0 0 0 0 0 0 1];
%TABLES

f=figure;
tabgp = uitabgroup(f,'Position',[.04 .5 .95 .45]);
cnames = {'dt=1/64','dt=1/128','dt=1/256','dt=1/512', 'dt=1/1024', 'dt=1/2048', 'dt=1/4096'};
rnames = {'N=3', 'N=7','N=15', 'N=31'};

tab = uitab(tabgp,'Title','Stable Solutions');
data = [stable];
% Create the uitable
t = uitable(tab,'Data', data,'ColumnName', cnames,'RowName', rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);




% Worksheet 4
% (E. Breznik, J. C. Medina, I. Tominec)
%
% Solving 2D unstationary PDE with space (FDM) discretization and Explicit/Implicit Euler method.
% -----------------------------------------------
clear all; close all;

% Setting the Simulation parameters. Grid size, time step, initial values.
NxNy=[3, 7, 15, 31]; % Grid size Nx=Ny for different cases.
dt=[1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 1/4096]; % Delta t's for different cases.
times=[1/8, 2/8, 3/8, 4/8]; % At these times we want to observe the solution of PDE.
y0=1; % initial value u(x,y,0)=1. It is the same for each value M(i,j) in the grid.

tic;
for j=1:length(times)
	
	folder=['t=' int2str(j) ':8']; % Making a folder for images within the current time at which we observe PDE.
	mkdir(folder);

	figure; title(['The solution plot for T(x,y,1/' int2str(times(j)) ').']); % Let's make the main windowz for plots.
	plotCnt=1; % The counter for subplots in one window.

	for i=1:length(NxNy) % Calculation for different NxNy.
		hx=1/(NxNy(i)+1); % Discretization in x direction.
		hy=hx; % Discretization in y direction.

		for k=1:length(dt) % Calc. for different delta t's.
			
			 % Preparing the matrix for solutions.
			M = zeros(NxNy(i), NxNy(i));
                        M(2:end-1, 2:end-1)=y0; % The initial value T(x,y,0)=y0 inside, and Dirichlet on the fense.

			 %% Solving
			 % ! Since the Grid is also time dependent, we have to change the RHS for each t+dt(k), until we arrive to the final value of u(x,y,times(j)).
			for o=1:times(j)/dt(k)+1
				for l=2:length(M)-1
					for m=2:length(M)-1
							rhs=(M(l+1,m)-2*M(l,m)+M(l-1,m))/hx^2 + (M(l,m+1)-2*M(l,m)+M(l,m-1))/hy^2;
							M(l,m)=Euler(rhs, dt(k), M(l,m), times(j)); % We just make one Explicit Euler step at a time. Next initial value y0 will be the value of M_{l, m-1}
					end
				end
			end
			
			 % Plotting (we are inside the k for loop).			
			subplot(length(NxNy), length(dt), plotCnt); % For each dt and current NxNy setting we make a subplot.
                        plotCnt=plotCnt+1;
			mesh(M); % We plot the solution to the set subplot area.

			figure(99); mesh(M); % We plot the solution to separate figure, which we also save as a PNG to folder and file with meaningful names.
			imgPath=[folder '/' int2str(NxNy(i)) 'x' int2str(NxNy(i)) '-1:' int2str(2^(5+k)) '.png'];
			title(imgPath); saveas(99, imgPath);
			close 99
		end
	end
end
toc
