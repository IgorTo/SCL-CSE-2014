function res=AdamsL1(dt, y0, tend)


res=zeros(1, tend/dt);
res(1)=y0; % Initial value.

P=@(y0, dt)(y0+dt/2*(14*y0-7/10*y0^2))/((y0*7*dt/20)+1); % The linearisation #1. We do not need the Newton method.

for i=1:length(res)-1
    res(i+1)=P(res(i), dt);
end
 
end
