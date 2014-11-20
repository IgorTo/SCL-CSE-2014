function [res, t]=Heun(func, dt, y0, tend)
% Heun method for solving ODE's.
%       Heun(f, s, y0, d) returns values of the
%       solution function (for ode y'=f(y,t)) at points 0:s:d.
%       f ... anonymous function, left side of diff.equation
%       y0 ... initial value

t=0:dt:tend;

res=zeros(1, length(t));
res(1)=y0;

for i=1:length(res)-1
	res(i+1)=res(i)+dt/2*(func(res(i), t(i))+func(res(i)+dt*func(res(i), t(i)), t(i)+dt));
end

end
