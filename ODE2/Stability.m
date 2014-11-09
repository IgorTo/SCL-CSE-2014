% STABILITY(values, time) determines whether the given case of simulation is stable.
%	Criteria:
%	1. The method for the given problem can be stable if the array of solutions is of the same length as the array of time at the given dt. Our simulation usually stops when certain values can not be found. This means that length(values) < length(time). We claim this happening as a contribution to non-stability of the method.
%	2. The method can be stable if the given data array appears to be convergent. We check that by criteria of abs(x_n - x_{n-1}) < eps, where n is the index of the final calculated value. For our problem, we used the eps of 1e-1.
%	If 1. AND 2. is true, then the simulation within certain method is stable.
%
%	Usage: Stability(values, time) 
%	values ... the data array
%	time ... the time array 0:dt:tend
%
%	The function returns '0' for non-stability and '1' for stability.
function res=Stability(values, time)

if (length(values) == length(time))
	if(abs(values(end)-values(end-1)) < 1e-1)
		res=1;
	else 
		res=0; 
	end
else
	res=0;
end
end
