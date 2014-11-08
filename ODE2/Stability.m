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
