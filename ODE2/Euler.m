function [res, t]=Euler(func, dt, y0, tend)
% Euler method for solving ODE's.
% 	Euler(f, dt, y0, tend) returns values of the
% 	solution function (for ode y'=f(y,t)) at points 0:dt:tend.
%	f ... anonymous function, left side of diff.equation
% 	y0 ... initial value

t=0:dt:tend; 
res=zeros(1, length(t));
res(1)=y0;

for i=1:length(res)-1
	res(i+1)=res(i)+func(res(i), t(i))*dt;
end

end
