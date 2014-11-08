function [res,t]=AdamsL2(func, dt, y0, tend)
t=0:dt:tend;

res=zeros(1, length(t));
res(1)=y0;

x0=1;
G=@(p,x0,dt) x0-p-dt/2*(7*p-7*p^2/10+7*x0-7*x0*p/10);
dG=@(p,x0,dt) 1-dt/2*(7-7*p/10);
for i=1:length(res)-1
    res(i+1)=res(i)+dt*.5*(func(res(i) ,t)+ 7*(1-res(i)/10)*Newton(G,dG,res(i),x0,dt)); 
   
end

end

function dp= Pprime(p,t)    %Calculates de derivative for a certain p and time
   dp=7*(1-p./10).*p;
end

function p = f_ex(t)        %Exact solution of the ODE at a certain time
    p=200./(20-10.*exp(-7.*t));
end