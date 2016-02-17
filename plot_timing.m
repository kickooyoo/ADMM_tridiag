colors = 'cmbrgky';
markers = 'ox+sd*.^v><ph';
exps = who('time*');
start_ndx = 1;
end_ndx = niters+1;
xform = @(x) 20*log10(x);
figure; hold on; 
legend_str = {};
for ii = 1:length(exps)
	curr_name = exps{ii};
	suffix_ndx = strfind(curr_name, '_');
	suffix = curr_name(suffix_ndx + 1 : end);
	if eval(sprintf('prod(size(%s)) ~= niters +1', curr_name))
		Nd = eval(sprintf('size(%s, 2)', curr_name));
		colm = true;
	else
		Nd = 1;
		colm = false;
	end
	curr_color = colors(mod(ii, length(colors)) + 1);
	for jj = 1:Nd
		curr_marker = markers(mod(jj, length(markers)) + 1);
		if colm
			eval(sprintf('plot(cumsum(%s(start_ndx:end_ndx,jj)), xform(nrmsd_%s(start_ndx:end_ndx,jj)), ''%s%s'')', curr_name, suffix, curr_color, curr_marker));
		else
			eval(sprintf('plot(cumsum(%s(:,start_ndx:end_ndx)), xform(nrmsd_%s(:,start_ndx:end_ndx)), ''%s'')', curr_name, suffix, curr_color));
		end
		hold on; 
		legend_str = [legend_str; sprintf('%s,%d', suffix,jj)];
	end
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
	eval(sprintf('plot(xform(nrmsd_%s(start_ndx:end_ndx)), ''%s'')', suffix, colors(mod(ii, length(colors)) + 1)));
	hold on; 
	legend_str = [legend_str; suffix];
end
legend(legend_str)
xlabel('iteration number')
title(sprintf('NRMSD to x^* over iteration for [%d %d]', Nx, Ny))
ylabel('NRMSD to x^* (dB)');
