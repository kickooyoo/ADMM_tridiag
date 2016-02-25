colors = 'cmgkbry';
markers = 'ox+sd*.^v><ph';
exps = who('time*');
start_ndx = 1;
end_ndx = niters+1;
xform = @(x) 20*log10(x);

if ~isvar('FontSize')
	FontSize = 12;
end

y_val = 'nrmsd'; %'nrmsd'; %'costOrig'

fhandle = figure('Position', [100 100 400 300]); hold on; 
legend_str = {};
for ii = 1:length(exps)
	curr_name = exps{ii};
	suffix_ndx = strfind(curr_name, '_');
	suffix = curr_name(suffix_ndx : end);
	if eval(sprintf('prod(size(%s)) ~= niters +1', curr_name))
		Nd = eval(sprintf('prod(size(%s))/(niters + 1)', curr_name));
	else
		Nd = 1;
	end
	if eval(sprintf('size(%s, 1) == niters + 1', curr_name))
		indeces = 'start_ndx:end_ndx,jj';
	else
		indeces = 'jj,start_ndx:end_ndx';
	end
	curr_color = colors(mod(ii, length(colors)) + 1);
	for jj = 1:Nd
		curr_marker = markers(mod(jj, length(markers)) + 1);
		eval(sprintf('plot(cumsum(%s(%s)), xform(%s%s(%s)), ''%s%s'')', curr_name, indeces, y_val, suffix, indeces, curr_color, curr_marker));
		legend_str = [legend_str; sprintf('%s,%d', suffix,jj)];
		hold on; 
	end
end
legend(legend_str, 'FontSize', FontSize)
xlabel('wall time (s)', 'FontSize', FontSize)
title(sprintf('%s to x^* over time for [%d %d]', upper(y_val), Nx, Ny))
ylabel('NRMSD to x^* (dB)', 'FontSize', FontSize);
axis tight

fhandle = figure('Position', [100 100 400 300]); hold on; 
legend_str = {};
for ii = 1:length(exps)
	curr_name = exps{ii};
	suffix_ndx = strfind(curr_name, '_');
	suffix = curr_name(suffix_ndx : end);
	if eval(sprintf('prod(size(%s)) ~= niters +1', curr_name))
		Nd = eval(sprintf('prod(size(%s))/(niters + 1)', curr_name));
	else
		Nd = 1;
	end
	if eval(sprintf('size(%s, 1) == niters + 1', curr_name))
		indeces = 'start_ndx:end_ndx,jj';
	else
		indeces = 'jj,start_ndx:end_ndx';
	end
	curr_color = colors(mod(ii, length(colors)) + 1);
	for jj = 1:Nd
		curr_marker = markers(mod(jj, length(markers)) + 1);
		eval(sprintf('plot(xform(%s%s(%s)), ''%s%s'')', y_val, suffix, indeces, curr_color, curr_marker));
		legend_str = [legend_str; sprintf('%s,%d', suffix,jj)];
		hold on; 
        end
end
legend(legend_str, 'FontSize', FontSize)
xlabel('iteration number', 'FontSize', FontSize)
title(sprintf('%s to x^* over iteration for [%d %d]', upper(y_val), Nx, Ny))
ylabel('NRMSD to x^* (dB)', 'FontSize', FontSize);
