% ADAMSL1(dt, y0, tend) solves the ODE equation of p'(t)=(1-p(t)/10)*7p(t)
%	The solver produces the results via the first linearised version of Adams-Moulton's method.
%	Usage: AdamsL1(dt, y0, tend)
%	dt ... time step
%	y0 ... initial value at time=0
%	tend ... the end time of simulation

function res=AdamsL1(dt, y0, tend)

res=zeros(1, tend/dt+1);
res(1)=y0; % Initial value.

P=@(y0, dt)(y0+dt/2*(14*y0-7/10*y0^2))/((y0*7*dt/20)+1); % The linearisation #1. We do not need the Newton method.

for i=1:length(res)-1
    res(i+1)=P(res(i), dt);
end
 
end
