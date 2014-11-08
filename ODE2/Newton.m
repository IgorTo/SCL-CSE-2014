function zero = Newton(f, df, t, y0)
% NEWTON(f, df, y0, dt) uses Newtons method 
% to compute a zero of a function f(p, t, y0)
% with respect to p.
%   df ... derivative of f
%   t .... time, at which initial approx. was computed
%   y0 ... initial approximation value

yold = y0+1;
ynew = y0; %starting approx. for the new value we are trying to compute
ok = 1;
iterat = 0;

while abs(yold-ynew)>1e-4
    yold = ynew;

    %stopping criteria (if eq. is not solvable):
    if (df(yold, t, y0) == 0 || iterat > 2000)
        ok = 0;
	disp(df(yold, t, y0))
	disp(iterat)
        break;
    else
        ynew = yold - f(yold, t, y0)/df(yold, t, y0);
    end
    iterat=iterat+1;
end

if ok==0
    zero = NaN;
else
    zero = ynew;
end

end
