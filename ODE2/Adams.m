function res = Adams(f, df, y0, dt, tend)
%ADAMS(f,dt,y0,tend) solves ODE y' = f(t, y)
%  with implicit Adams Moulton method.
%   f ... function (right side of ODE)
%  df ... derivative of f
%  y0 ... initial value
%  dt ... time step size
% tend .. end of observed time enterval

t = 0:dt:tend;

res = zeros(1, length(t));
res(1) = y0;

G = @(p,time,p0) p-p0-dt/2*(f(time, p0)+f(time+dt, p));
dG = @(p,time,p0) 1-dt/2*df(time+dt, p);
for i=1:length(res)-1
    res(i+1)=res(i)+dt*.5*(func(res(i) ,t)+func(Newton(G,dG,res(i),x0,dt),t+dt)); 
    
end
 
end

