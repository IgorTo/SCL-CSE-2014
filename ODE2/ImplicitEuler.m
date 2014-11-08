function  y = ImplicitEuler(f, df, y0, dt, tend)
% IMPLICITEULER(f, df, y0, dt, tend) uses implicit Euler
% method to solve ODE y' = f(t, y). Returns a vector of
% approximations of y(t).
%   f  .... function of y and t, left side of ODE
%   df .... derivative of f with respect to y
%   y0 .... initial value
%   dt .... time step
%  tend ... end of observed time interval


y = [y0]; %assumption: t0 = 0
steps = tend/dt;

G = @(p, time, start) (p-dt*f(time+dt, p)-start);
dG = @(p, time, start) (1-dt*df(time+dt, p));
for k=1:steps
    N = Newton(G, dG, dt*(k-1), y(k));
    %check for stopping criteria (if eq is not solvable):
    if isnan(N)
        break;
    else
        y(k+1) = N;
    end
    %disp(y);
end

end

    
    