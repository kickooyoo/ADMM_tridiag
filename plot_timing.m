colors = 'cmgkbry';
markers = 'ox+sd*.^v><ph';
exps = who('time*');
start_ndx = 1;
end_ndx = niters+1;
xform = @(x) 20*log10(x);

if ~isvar('FontSize')
	FontSize = 12;
end

fhandle = figure('Position', [100 100 400 300]); hold on; 
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
			eval(sprintf('plot(cumsum(%s(jj,start_ndx:end_ndx)), xform(nrmsd_%s(jj,start_ndx:end_ndx)), ''%s%s'')', curr_name, suffix, curr_color, curr_marker));
			legend_str = [legend_str; sprintf('%s,%d', suffix,jj)];
		else
			eval(sprintf('plot(cumsum(%s(:,start_ndx:end_ndx)), xform(nrmsd_%s(:,start_ndx:end_ndx)), ''%s'')', curr_name, suffix, curr_color));
			legend_str = [legend_str; sprintf('%s', suffix)];
		end
		hold on; 
	end
end
legend(legend_str, 'FontSize', FontSize)
xlabel('wall time (s)', 'FontSize', FontSize)
title(sprintf('NRMSD to x^* over time for [%d %d]', Nx, Ny))
ylabel('NRMSD to x^* (dB)', 'FontSize', FontSize);
axis tight

fhandle = figure('Position', [100 100 400 300]); hold on; 
legend_str = {};
for ii = 1:length(exps)
	curr_name = exps{ii};
	suffix_ndx = strfind(curr_name, '_');
	suffix = curr_name(suffix_ndx + 1 : end);
	eval(sprintf('plot(xform(nrmsd_%s(start_ndx:end_ndx)), ''%s'')', suffix, colors(mod(ii, length(colors)) + 1)));
	hold on; 
	legend_str = [legend_str; suffix];
end
legend(legend_str, 'FontSize', FontSize)
xlabel('iteration number', 'FontSize', FontSize)
title(sprintf('NRMSD to x^* over iteration for [%d %d]', Nx, Ny))
ylabel('NRMSD to x^* (dB)', 'FontSize', FontSize);
