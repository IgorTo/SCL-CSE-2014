close all; clear all;

% The problem definition.
p_der=@(p, t)(7*(1-p./10).*p);
p0=20;

% The Analytical solution.
p_ref=@(t)(200./(20-10.*exp(-7.*t)));

% A)
t=0:0.01:5; 
figure(1); plot(t, p_ref(t));

% b)
t0=0;
tend=5; % The Simulation parameters.
dt=1; % The starting delta t.

colors=['rbcygmk']; % Let's prepare the string of color codes for plotting.

% The Explicit Euler figure
figure(2); plot(t, p_ref(t), sprintf('%c',colors(1))); 
axis([t0 tend 0 20]); xlabel('time'); ylabel('population'); title('Explicit Euler');

% The Heun figure
figure(3); plot(t, p_ref(t), sprintf('%c',colors(1))); 
axis([t0 tend 0 20]); xlabel('time'); ylabel('population'); title('Heun');

% Solving the ODE and plotting its solution for different dt
for i=1:6
    [p, time]=Euler(p_der, dt, p0, tend);
    figure(2); hold on; plot(time, p, sprintf('%c', colors(i+1)));
    errorEuler(i)=sqrt(sum((p_ref(t0:dt:tend)-p).^2).*dt./5);
    
    [p, time]=Heun(p_der, dt, p0, tend);
    figure(3); hold on; plot(time, p, sprintf('%c',colors(i+1)));
    errorHeun(i)=sqrt(sum((p_ref(t0:dt:tend)-p).^2).*dt./5);
    dt=dt/2;
end

% The Plot Legend configuration
figure(2); legend('Analytic', 'dt=1', 'dt=1/2', 'dt=1/4', 'dt=1/8', 'dt=1/16', 'dt=1/32', 'Location', 'SouthEast');
figure(3); legend('Analytic', 'dt=1', 'dt=1/2', 'dt=1/4', 'dt=1/8', 'dt=1/16', 'dt=1/32', 'Location', 'SouthEast');
close all

%H)

%Reduced error
for i=1:5
    if(errorEuler(i)<exp(30))  %Taking care of divisions by infinity. 
         errorEuler_Red(i)=errorEuler(i)/errorEuler(i+1);
    else 
         errorEuler_Red(i)=NaN;   %If there is no reduced error to compute, a NaN appears on the table
    end
    
    if(errorHeun(i)<exp(30))
         errorHeun_Red(i)=errorHeun(i)/errorHeun(i+1);
    else 
         errorHeun_Red(i)=NaN;
    end
   
    %To include reduced errors for new methods
end



%Tables 


f=figure;
tabgp = uitabgroup(f,'Position',[.05 .6 .95 .3]);
cnames = {'dt=1/2','dt=1/4','dt=1/8','dt=1/16','dt=1/32'};
rnames = {'Error','Error red.'};

tab1 = uitab(tabgp,'Title','Euler implicit');

data = [errorEuler(2:6); errorEuler_Red(1:5)];
% Create the uitable
t = uitable(tab1,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);

tab2 = uitab(tabgp,'Title','Heun');
data = [errorHeun(2:6); errorHeun_Red(1:5)];
% Create the uitable
t = uitable(tab2,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);

%To include tables for new methods

