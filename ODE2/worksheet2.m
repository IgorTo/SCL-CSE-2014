close all; clear all;

% The problem definition.
p_der = @(p, t)(7*(1-p./10).*p);
dp_der = @(p, t)(7 - 7*p./5);
p0=20;

% The Analytical solution.
p_ref=@(t)(200./(20-10.*exp(-7.*t)));

% A)
t=0:0.01:5; 
figure(1); plot(t, p_ref(t));

% B)
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

% The Implicit Euler figure
figure(4); plot(t, p_ref(t), sprintf('%c',colors(1)));
axis([t0 tend 0 20]); xlabel('time'); ylabel('population'); title('ImplicitEuler');

% The Adams figure
figure(5); plot(t, p_ref(t), sprintf('%c',colors(1))); 
axis([t0 tend 0 20]); xlabel('time'); ylabel('population'); title('Adams');

% The AdamsL1 figure
figure(6); plot(t, p_ref(t), sprintf('%c',colors(1))); 
axis([t0 tend 0 20]); xlabel('time'); ylabel('population'); title('AdamsL1');

% The AdamsL2 figure
figure(7); plot(t, p_ref(t), sprintf('%c',colors(1))); 
axis([t0 tend 0 20]); xlabel('time'); ylabel('population'); title('AdamsL2');

% Initialization of error arrays
errorEuler = zeros(1,6);
errorHeun = zeros(1,6);
errorImpEuler = zeros(1,6);
errorAdams = zeros(1,6);
errorAdamsL1 = zeros(1,6);
errorAdamsL2 = zeros(1,6);

% Init of stability arrays
stabilityHeun=zeros(1,6);
stabilityEuler=zeros(1,6);
stabilityImpEuler=zeros(1,6);
stabilityAdams=zeros(1,6);
stabilityAdamsL1=zeros(1,6);
stabilityAdamsL2=zeros(1,6);

% Solving the ODE and plotting its solution for different dt
for i=1:6
    time_dt = t0:dt:tend;
    
    % Euler
    [p, time] = Euler(p_der, dt, p0, tend);
    figure(2); hold on; plot(time, p, sprintf('%c', colors(i+1)));
    errorEuler(i) = sqrt(sum((p_ref(time_dt)-p).^2).*dt./5);
    stabilityEuler(i)=Stability(p, time_dt);
    
    
    % Heun
    [p, time] = Heun(p_der, dt, p0, tend);
    figure(3); hold on; plot(time, p, sprintf('%c',colors(i+1)));
    errorHeun(i)=sqrt(sum((p_ref(time_dt)-p).^2).*dt./5);
    stabilityHeun(i)=Stability(p, time_dt);
    
    
    % Implicit Euler
    p = ImplicitEuler(p_der, dp_der, p0, dt, tend);
    success = length(p)-1;
    time = t0:dt:success*dt;
    figure(4); hold on; plot(time, p, sprintf('%c',colors(i+1)));
    errorImpEuler(i)=sqrt(sum((p_ref(time)-p).^2).*dt./5);
    stabilityImpEuler(i)=Stability(p, time_dt);
    
    
    % Adams
    p = Adams(p_der, dp_der, p0, dt, tend);
    success = length(p)-1;
    time = t0:dt:success*dt;
    figure(5); hold on; plot(time, p, sprintf('%c',colors(i+1)));
    errorAdams(i)=sqrt(sum((p_ref(time)-p).^2).*dt./5);
    stabilityAdams(i)=Stability(p, time_dt);
    
    
    % Adams Linearisation #1
    p = AdamsL1(dt, p0, tend); % Solution
    figure(6); hold on; plot(time_dt, p, sprintf('%c',colors(i+1))); % Plot
    errorAdamsL1(i)=sqrt(sum((p_ref(time_dt)-p).^2).*dt./5); % Error
    stabilityAdamsL1(i)=Stability(p, time_dt);
   
    % Adams Linearisation #2
    p = AdamsL2(dt, p0, tend); % Solution
    figure(7); hold on; plot(time_dt, p, sprintf('%c',colors(i+1))); % Plot
    errorAdamsL2(i)=sqrt(sum((p_ref(time_dt)-p).^2).*dt./5); % Error
    stabilityAdamsL2(i)=Stability(p, time_dt);
    
    
    % go to the next dt
    dt = dt/2;
end

% The Plot Legend configuration
figure(2); legend('Analytic', 'dt=1', 'dt=1/2', 'dt=1/4', 'dt=1/8', 'dt=1/16', 'dt=1/32', 'Location', 'SouthEast');
figure(3); legend('Analytic', 'dt=1', 'dt=1/2', 'dt=1/4', 'dt=1/8', 'dt=1/16', 'dt=1/32', 'Location', 'SouthEast');
figure(4); legend('Analytic', 'dt=1', 'dt=1/2', 'dt=1/4', 'dt=1/8', 'dt=1/16', 'dt=1/32', 'Location', 'SouthEast');
figure(5); legend('Analytic', 'dt=1', 'dt=1/2', 'dt=1/4', 'dt=1/8', 'dt=1/16', 'dt=1/32', 'Location', 'SouthEast');
figure(6); legend('Analytic', 'dt=1', 'dt=1/2', 'dt=1/4', 'dt=1/8', 'dt=1/16', 'dt=1/32', 'Location', 'SouthEast');
figure(7); legend('Analytic', 'dt=1', 'dt=1/2', 'dt=1/4', 'dt=1/8', 'dt=1/16', 'dt=1/32', 'Location', 'SouthEast');


%H)
%Reduced error
errorEuler_Red = zeros(1,6);
errorHeun_Red = zeros(1,6);
errorImpEuler_Red = zeros(1,6);
errorAdams_Red = zeros(1,6);
errorAdamsL1_Red = zeros(1,6);
errorAdamsL2_Red = zeros(1,6);
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
   
    errorImpEuler_Red(i) = errorImpEuler(i)/errorImpEuler(i+1);
    errorAdams_Red(i)=errorAdams(i)/errorAdams(i+1);
    errorAdamsL1_Red(i)=errorAdamsL1(i)/errorAdamsL1(i+1);
    errorAdamsL2_Red(i)=errorAdamsL2(i)/errorAdamsL2(i+1);
end



%Tables 


f=figure;
tabgp = uitabgroup(f,'Position',[.05 .6 .95 .3]);
cnames = {'dt=1/2','dt=1/4','dt=1/8','dt=1/16','dt=1/32'};
rnames = {'Error','Error red.'};

tab1 = uitab(tabgp,'Title','Explicit Euler');
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

tab3 = uitab(tabgp,'Title','Implicit Euler');
data = [errorImpEuler(2:6); errorImpEuler_Red(1:5)];
% Create the uitable
t = uitable(tab3,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);

tab4 = uitab(tabgp,'Title','Adams');
data = [errorAdams(2:6); errorAdams_Red(1:5)];
% Create the uitable
t = uitable(tab4,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);

tab5 = uitab(tabgp,'Title','AdamsL1');
data = [errorAdamsL1(2:6); errorAdamsL1_Red(1:5)];
% Create the uitable
t = uitable(tab5,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);

tab6 = uitab(tabgp,'Title','AdamsL2');
data = [errorAdamsL2(2:6); errorAdamsL2_Red(1:5)];
% Create the uitable
t = uitable(tab6,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);

% Stable cases Plot

cnames = {'dt=1/2','dt=1/4','dt=1/8','dt=1/16','dt=1/32'};
rnames = {'Explicit Euler', 'Heun', 'Implicit Euler', 'Adams-Moulton', 'Adams-Moulton L1', 'Adams-Moulton L2'};
tab7 = uitab(tabgp,'Title','Stable cases');
data = [stabilityEuler(2:6); stabilityHeun(2:6); stabilityImpEuler(2:6); stabilityAdams(2:6); stabilityAdamsL1(2:6); stabilityAdamsL2(2:6)];
% Create the uitable
t = uitable(tab7,'Data',data,'ColumnName',cnames,'RowName',rnames);
% Set width and height
t.Position(3) = t.Extent(3);
t.Position(4) = t.Extent(4);


