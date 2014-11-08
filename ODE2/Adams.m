function res = Adams(f, df, y0, dt, tend)
%ADAMS(f,dt,y0,tend) solves ODE y' = f(y, t)
%  with implicit Adams Moulton method.
%   f ... function (right side of ODE)
%  df ... derivative of f
%  y0 ... initial value
%  dt ... time step size
% tend .. end of observed time enterval

t = 0:dt:tend;
l = length(t);
res = zeros(1, l); %l=steps+1
res(1) = y0;

G = @(p,time,p0) (p-p0-dt/2*(f(p0, time)+f(p, time+dt)));
dG = @(p,time,p0) (1-dt/2*df(p, time+dt));

for i=1:l-1
    N = Newton(G, dG, t(i), res(i));
	disp(N)
    if isnan(N)
        break;
    else
        res(i+1) = N;
    end
end
 
end

