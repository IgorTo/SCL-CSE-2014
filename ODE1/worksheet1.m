%                  WORKSHEET 1- EXPLICIT ODE (E. Breznik, J. C. Medina, I. Tominec)
%
% This program solves an ODE for the following explicit methods:
% explicit Euler, Heun and Runge-Kutta. The ODE is p'=(1-p/10)*p
% Moreover, the change of the error by reducing the time steps 
% for 1,1/2,1/4 and 1/8 and also an approximation error are calculated 
%
% Note:
% Design considerations: *The graphics for each dt are shown sepparately in
% the same subplot to appreciate better the change of accuracy. 
% *Each method must be chosen in the code by uncommenting it (@ line 37). We chose not to
% run the three methods simultaniously so that there would be no 6
% different graph/table windows at the same time.
%                         26 October 2014


function worksheet1

 %Initial conditions
 t_start=0;
 t_end=5;
 initval=1;

 DT=[1/8,1/4,1/2,1,1/16];  % We need  the 1/16 value for de reduced error
 
 error1=zeros(1,5);      
 errorRed=zeros(1,4);
 errorApp=zeros(1,4);
 
 
 
 for i=1:5
  f=[initval];           % We begin the function with the initial condition
  dt=DT(i);        % A different dt is chosen for each for cycle 
  N=(t_end-t_start)/dt;   
  
  
% !!! MANUALLY please choose between ExplicitEuler, Heun or RungeKutta 
  [f,error]=ExplicitEuler(f,t_start,dt,N);
  %[f,error]=Heun(f,t_start,dt,N);
  %[f,error]=RungeKutta(f,t_start,dt,N);
  
  
  error1(i)=error;
        
 if(i==5)    % This is the case when dt=1/16 and there is no neccesity for plotting it
     break
 end
 
 %%%%%%%APPROXIMATE ERROR%%%%%%%%%%%
 switch i
     case 1 
        g=f;        % G is the array that stores the values of the best approximation
     case 2
         for j=2:N+1                                   % j=1 is the initial conditions. It runs until N+1 because of the N steps plus the initial condition.  
         errorApp(2)=errorApp(2)+(f(j)-g(2*j-1))^2;    % Each value of the dt=1/4  is compared for the respective 2*j-1 of g
         end
        errorApp(2)=sqrt(errorApp(2)*dt/5);
     case 3
         for j=2:N+1
          errorApp(3)=errorApp(3)+(f(j)-g(4*j-3))^2;  % Each value of the dt=1/2  is compared for the respective 4*j-2 of g
         end
         errorApp(3)=sqrt(errorApp(3)*dt/5); 
     case 4
         for j=2:N+1
         errorApp(4)= errorApp(4)+(f(j)-g(8*j-7))^2;   % Each value of the dt=1 is compared for the respective 2*j-1 of g
         end
          errorApp(4)=sqrt( errorApp(4)*dt/5);    
       
 end   

%%%%%%PLOTTING%%%%%%%%%%%%% 
 
 % Plot numerical solution (red)
subplot(2,2,i)
 plot(t_start + dt*[0:N], f, 'r*');
hold on

% Compare with exact solution (blue)
t_plot = [ t_start : 0.01*dt : t_end ];
plot(t_plot, f_ex(t_plot), 'b');

%Plot design
xlabel('time');
ylabel('population');
title({'dt=',dt});
xlabel('Time');

clear f    % Beginn f again from scratch for next dt
 end       % END of the four iterations

%%%%%%%REDUCED ERROR%%%%%%%
errorRed(1)=error1(1)/error1(5);  % We need the extra value for dt=1/16
errorRed(2)=error1(2)/error1(1);
errorRed(3)=error1(3)/error1(2);
errorRed(4)=error1(4)/error1(3);
 
 
 
%%%%%%%%CREATE THE TABLE%%%%%%%%%%
 f = figure('Position',[440 500 461 146]);

% create the data
data = [error1(1:4); errorRed(1:4);errorApp(1:4)];

% Create the column and row names in cell arrays 
cnames = {'dt=1/8','dt=1/4','dt=1/2','dt=1'};
rnames = {'Error','Error red.','Error app.'};

% Create the uitable
t = uitable(f,'Data',data,...
            'ColumnName',cnames,... 
            'RowName',rnames);

% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);

end

%%%%NUMMERICAL METHODS AS FUNCTIONS%%%%

function [f,error]=ExplicitEuler(f,t_start,dt,N)
error=0;
for t = t_start + dt*( 1:(N) )
    f(end+1)=f(end)+dt*Pprime(f(end) ,t);   % Instead of using an index k to save the values of Yn+1, we use the f(end+1)=F(f(end)) to create the array 
    error=error+(f(end)-f_ex(t))^2;
    
end
error=sqrt(error*dt/5);
end

function [f,error]=Heun(f,t_start,dt,N)
error=0;
for t = t_start + dt*( 1:(N) )
    f(end+1)=f(end)+dt*.5*(Pprime(f(end) ,t)+Pprime(f(end)+dt*Pprime(f(end),t),t+dt)); 
    error=error+(f(end)-f_ex(t))^2;
end
 error=sqrt(error*dt/5);
end

function [f,error]=RungeKutta(f,t_start,dt,N)
error=0;
for t = t_start + dt*( 1:(N) )
    P1=Pprime(f(end),t);
    P2=Pprime(f(end)+(dt/2)*P1,t+dt/2);
    P3=Pprime(f(end)+(dt/2)*P2,t+dt/2);
    P4=Pprime(f(end)+dt*P3,t+dt);
    f(end+1)=f(end)+(dt/6)*(P1+2*P2+2*P3+P4);
    error=error+(f(end)-f_ex(t))^2;          
end
error=sqrt(error*dt/5);
end

%%%%DERIVATIVE AND EXACT FUNCTION%%%%%
function dp= Pprime(p,t)    % Calculates de derivative for a certain p and time
   dp=(1-p/10)*p;
end

function p = f_ex(t)        % Exact solution of the ODE at a certain time
    p=10./(1+9*exp(-t));
end







