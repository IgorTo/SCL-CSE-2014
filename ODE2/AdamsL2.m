% ADAMSL2(dt, y0, tend) solves the ODE equation of p'(t)=(1-p(t)/10)*7p(t)
%       The solver produces the results via the second linearised version of Adams-Moulto$
%       Usage: AdamsL2(dt, y0, tend)
%       dt ... time step
%       y0 ... initial value at time=0
%       tend ... the end time of simulation

function res=AdamsL2(dt, y0, tend)

res=zeros(1, tend/dt+1);
res(1)=y0; % Initial value.

P=@(y0, dt)(y0+dt/2*(7*y0-7/10*y0^2))/(1+dt/2*(7*y0/10-7)); % The linearisation #2. We do not need the Newton method -> the polynomial's order is 1.

for i=1:length(res)-1
    res(i+1)=P(res(i), dt);
end
end
