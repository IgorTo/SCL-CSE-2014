function res=Euler(rhs, dt, y0, tend)
% Euler method for solving ODE's.
% 	Euler(f, dt, y0, tend) returns values of the
% 	solution function (for ode y'=f(y,t)) at points 0:dt:tend.
%	rhs ... anonymous function, left side of diff.equation
% 	y0 ... initial value
%	Note: The function is implemented to make just one Euler method step.
%	      Although it may work for other kinds of problems, this function
%	      is adjusted specifically for solving ODE derived from unstatio-
%             nary R^2 PDE heat equation.


res=y0+rhs*dt;

end
