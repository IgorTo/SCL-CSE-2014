%   Worksheet 4
% (E. Breznik, J. C. Medina, I. Tominec)
%
% -Solving time dependant partial differential equations.
% -----------------------------------------------
clear all; close all;

%a)Determine the function when the limit tends to zero

% We assumme that the temperature doesn't change and
% therefore the derivative with respect to time is equal to zero.
% We have the 2D Laplace equation and since there is no additional 
% heat source, the solution is zero everywhere.

%b) Implement an explicit Euler step 
% Function call to explicitEuler2D(Nx,Ny,dt,sol)

%c) Solve for the values of Nx,Ny and dt for times= 1/8,2/8,3/8 and 4/8
N = [3,7, 15, 31];


%plotting:
figure(1);
set(gcf,'numbertitle','off','name','Time=1/8');
figure(2);
set(gcf,'numbertitle','off','name','Time=2/8');
figure(3);
set(gcf,'numbertitle','off','name','Time=3/8');
figure(4);
set(gcf,'numbertitle','off','name','Time=4/8');
counter=1;


for i=1:4
    Nx = N(i);
    Ny = N(i);
    dT=1/64; %First dT to obtain.
    %this we need later for plotting:
    [x,y] = meshgrid(1/(Nx+1)*(0:(Nx+1)), 1/(Ny+1)*(0:(Ny+1)));
    
    f= ones(Nx*Ny,1); %Initial condition vector
    for deltaTimes=1:7
        sol1=f; %Initialize solution to time zero, for every dT
        repetitions=2^(deltaTimes-1)*8; %Depending on the dT, the number of iterations of the Euler differ 
        
        for times=1:4 
            
            %We use a matrix to make the program more efficient (Instead of
            %solving the ODE point by point)
            
           
            M=EulerMatrix(Nx,Ny,dT);
           for j=1:repetitions
                sol1=explicitEuler2D(Nx,Ny,dT,sol1,M);
           end
           
           figure(times);
           Temp = fliplr(reshape(sol1, [Ny,Nx]))';
           Temp = [zeros(1,Nx+2);zeros(Ny,1),Temp,zeros(Ny,1);zeros(1,Nx+2)];
           subplot(4,7,counter);
           surf(x,y,Temp);
           title({dT,Nx});
        end
        counter=counter+1;
        dT=dT/2;
       
    end
 
end

    
    
    