function beta = choose_beta(orient, slice, reduction)
%function beta = choose_beta(orient, slice, reduction)

switch orient
case 'sim'
	beta = choose_beta_sim(reduction);
case 'axial'
	beta = choose_beta_axial(slice, reduction);
case 'sagittal'
	beta = 2^23;
case 'coronal'
	beta = 2^21;
otherwise
	display(sprintf('unknown orientation %s', orient));
	keyboard;
end

end


function beta = choose_beta_sim(reduction)
switch reduction
case {4, 6}
	beta = 2^12;% cannot find file, need to regen
	beta = 2^13;
case 10
	beta = 2^14;
otherwise 
	display('unknown which is best beta, choosing 2^13 arbitrarily');
	beta = 2^13;
end
end

function beta = choose_beta_axial(slice, reduction)
switch slice
case 67
	beta = 2^27;
case 38
	if reduction <= 8
		beta = 2^20;
	elseif reduction == 12
		beta = 2^18;
	else
		beta = 2^20;
	end
case 90
	beta = 2^25;
otherwise
	display('unknown which is best beta, choosing 2^20 arbitrarily');
	beta = 2^20;
end
end
