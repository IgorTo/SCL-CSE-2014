function x1 = Newton( f, df, x0, eps )
%NEWTON A function for solving the equation F(x)=0 by Newton's method.
%   Usage: Newton(f, df, x0)
%   f ... left hand-side function F(x), should be given by anonymous
%   function
%   df ... the derivative of left hand-side function F'(x), should be given
%   by anonymous function
%   x0 ... initial (guess) value
%   eps

while(true)
    x1=x0-f(x0)/df(x0);
    if (abs(x1-x0) < eps) break; end
    x0=x1;
end
end
