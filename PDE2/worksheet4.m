% Worksheet 4
% (E. Breznik, J. C. Medina, I. Tominec)
%
% Solving 2D unstationary PDE with space (FDM) discretization and Explicit/Implicit Euler method.
% -----------------------------------------------
clear all; close all;


%a) Determine the function when the limit tends to zero

% We assumme that the temperature doesn't change in time and
% therefore the derivative with respect to time is equal to zero.
% We have the 2D Laplace equation and since there is no additional
% heat source, the solution is zero everywhere.

%% b+c) Explicit Euler

% Setting the Simulation parameters. Grid size, time step, initial values.
NxNy=[3, 7, 15, 31]; % Grid size Nx=Ny for different cases.
dt=[1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 1/4096]; % Delta t's for different cases.
times=[1/8, 2/8, 3/8, 4/8]; % At these times we want to observe the solution of PDE.

for j=1:length(times)

	% Saving + Plot manners.
	folder=['t_' int2str(j) ':8_']; % Making a folder for images within the current time at which we observe PDE.
	mkdir(folder);
	figure('Position', [100, 100, 1600, 900]); 
	axes('Position', [0, 0.95, 1, 0.05]);
	set(gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None');
	text(0.5, 0, ['Explicit Euler. The solution plot for T(x,y,t=' int2str(j) '/8).'], 'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');

	plotCnt=1; % The counter for subplots in one window.

	for i=1:length(NxNy) % Calculation for different NxNy.
		for k=1:length(dt) % Calc. for different delta t's.

			solution=ones(NxNy(i)*NxNy(i), 1); % Vector of initial values u(x,y,0)=1.
			M=EulerMatrix(NxNy(i), NxNy(i), dt(k)); % Tridiagonal matrix with Explicit Euler implemented.

			 %% Solving
			 % ! Since the Grid is also time dependent, we have to change the RHS for each t+dt(k), until we arrive to the final value of u(x,y,times(j)).
			for o=1:times(j)/dt(k)+1
				solution=explicitEuler2D(NxNy(i), dt(k), solution, M);
			end

			solution = fliplr(reshape(solution, [NxNy(i), NxNy(i)]))'; % Reshaping the vector into matrix of NxNy*NxNy.
		        solution = [zeros(1, NxNy(i)+2); zeros(NxNy(i),1), solution, zeros(NxNy(i),1); zeros(1,NxNy(i)+2)]; % Adding the Dirichlet on the fense.

			 % (Sub)Plotting.
			subplot(length(NxNy), length(dt), plotCnt); % For each dt and current NxNy setting we make a subplot.
                        plotCnt=plotCnt+1;
			mesh(solution); % We plot the solution to the set subplot area.
			title([int2str(NxNy(i)) 'x' int2str(NxNy(i)) ', dt=1/' int2str(2^(5+k))]);
			set(gca, 'FontSize', 9);

			 % Plotting each solution in separate window for storing purposes.
			figure(99); mesh(solution);
			imgPath=[folder '/' int2str(NxNy(i)) 'x' int2str(NxNy(i)) '-1:' int2str(2^(5+k)) '.png'];
			title([int2str(NxNy(i)) 'x' int2str(NxNy(i)) ', dt=1/' int2str(2^(5+k))]); saveas(99, imgPath);
			close 99
		end
	end
	saveas(j, [folder '/0SubPlots.png']); % We store the figure with many subplots.
end


%%% d+e) Implicit Euler
dt=1/64;

% Plot preparation.
f=figure('Position', [100, 100, 1600, 900]);
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, ['The plot for solutions made with Implicit Euler.'], 'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

% Calculation.
plotCnt=1;
for j=1:length(times)
	for k=1:length(NxNy)

		 % We prepare the solution Matrix with Dirichlet on the fense and u(x,y,0)=1 everywhere else.
		solution=zeros(NxNy(k)+2, NxNy(k)+2); solution(2:end-1, 2:end-1)=1; 

		 % Solving.
		for o=1:times(j)/dt+1
			solution=GausSeidel(NxNy(k), NxNy(k), solution, dt);
		end

		% Plotting.
		subplot(length(NxNy), length(times), plotCnt);
		plotCnt=plotCnt+1;
		mesh(solution)
		title([int2str(NxNy(k)) 'x' int2str(NxNy(k)) '@ t=1/' int2str(j)]);
		set(gca, 'FontSize', 9);
	end
end
saveas(f, 'ImplicitEuler.png'); % Plot save.



%% Table of stable cases.
data = [1 1 1 1 1 1 1;
        0 0 1 1 1 1 1;
        0 0 0 0 1 1 1;
        0 0 0 0 0 0 1];
f=figure;
tabgp = uitabgroup(f,'Position',[.05 .5 .95 .45]);
cnames = {'dt=1/64','dt=1/128','dt=1/256','dt=1/512', 'dt=1/1024', 'dt=1/2048', 'dt=1/4096'};
rnames = {'Nx=Ny=3', 'Nx=Ny=7','Nx=Ny=15', 'Nx=Ny=31'};
tab7 = uitab(tabgp,'Title','Stable cases @ Implicit Euler');
% Create the uitable
t = uitable(tab7,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);
