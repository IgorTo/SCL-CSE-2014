function x1 = Newton( f,df, p,x0, dt )
%NEWTON A function for solving the equation F(x)=0 by Newton's method.
%   Usage: Newton(f, df, x0, eps)
%   f ... left hand-side function F(x), should be given by anonymous
%   function
%   df ... the derivative of left hand-side function F'(x), should be given
%   by anonymous function
%   x0 ... initial (guess) value
%   eps


    n=0;
while(true)
    fprime=df(p,x0,dt);
    if(fprime==0|| n>100)
        fprintf('error');
        break;
    end
    x1=x0-f(p,x0,dt)/fprime;
    if (abs(x1-x0) < 1e-4) break; end
    x0=x1;
    n=n+1;
end
end
