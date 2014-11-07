function [res,t]=AdamsL1(func, dt, y0, tend)


t=0:dt:tend;

res=zeros(1, length(t));
res(1)=y0;

x0=1;
G=@(p,x0,dt) x0-p-dt/2*(7*p-7*p^2/10+7*p-7*x0*p/10);
dG=@(p,x0,dt) 1+dt/2*(7*p/10);
for i=1:length(res)-1
    res(i+1)=res(i)+dt*.5*(func(res(i) ,t)+ 7*(1-Newton(G,dG,res(i),x0,dt)/10)*res(i)); 
    
end
 
end
