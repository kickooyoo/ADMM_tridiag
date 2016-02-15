colors = 'cmbrgk';
exps = who('time*');
start_ndx = 2;
end_ndx = niters+1;
xform = @(x) 20*log10(x);
figure; hold on; 
legend_str = {};
for ii = 1:length(exps)
	curr_name = exps{ii};
	suffix_ndx = strfind(curr_name, '_');
	suffix = curr_name(suffix_ndx + 1 : end);
	eval(sprintf('plot(cumsum(%s(start_ndx:end_ndx)), xform(nrmsd_%s(start_ndx:end_ndx)), ''%s'')', curr_name, suffix, colors(ii)));
	hold on; 
	legend_str = [legend_str; suffix];
end
legend(legend_str)
xlabel('wall time (s)')
title(sprintf('NRMSD to x^* over time for [%d %d]', Nx, Ny))
ylabel('NRMSD to x^* (dB)');
axis tight

figure; hold on; 
legend_str = {};
for ii = 1:length(exps)
	curr_name = exps{ii};
	suffix_ndx = strfind(curr_name, '_');
	suffix = curr_name(suffix_ndx + 1 : end);
	eval(sprintf('plot(xform(nrmsd_%s(start_ndx:end_ndx)), ''%s'')', suffix, colors(ii)));
	hold on; 
	legend_str = [legend_str; suffix];
end
legend(legend_str)
xlabel('iteration number')
title(sprintf('NRMSD to x^* over iteration for [%d %d]', Nx, Ny))
ylabel('NRMSD to x^* (dB)');
