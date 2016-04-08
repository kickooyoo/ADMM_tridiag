colors = 'cmgkbry';
markers = 'ox+sd*.^v><ph';
<<<<<<< HEAD
% markers = '-';
=======
%markers = '-';
>>>>>>> 30b803c31ccfabd8a6f906123705ca3e682028a9
exps = who('time*');
start_ndx = 1;
end_ndx = niters+1;
if ~isvar('orn_ndx')
        orn_ndx = round(end_ndx/2);
end
xform = @(x) 20*log10(x);

if ~isvar('FontSize')
	FontSize = 12;
end

y_val = 'nrmsd'; %'nrmsd'; %'costOrig'
xlabels = {'wall time (s)'; 'iteration number'};
for plot_ndx = 1:2
fhandle = figure('Position', [100 100 400 300]); hold on; 
legend_str = {};
for ii = 1:length(exps)
	curr_name = exps{ii};
	suffix_ndx = strfind(curr_name, '_');
	suffix = curr_name(suffix_ndx : end);
        if eval(sprintf('size(%s,1) == niters + 1', curr_name))
                eval(sprintf('%s = %s.'';', curr_name, curr_name));
                eval(sprintf('%s%s = %s%s.'';', y_val, suffix, y_val, suffix));
        end
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
        if start_ndx == 1
                tcomp = sprintf(' - %s(1)', curr_name);
        end
	curr_color = colors(mod(ii, length(colors)) + 1);
	for jj = 1:Nd
                if plot_ndx == 1
                        plot_command = sprintf('h(ii,1) = plot(cumsum(%s(%s))%s, xform(%s%s(%s)), ''%s%s'');', curr_name, indeces, tcomp, y_val, suffix, indeces, curr_color, '-');
                else
                        plot_command = sprintf('h(ii,1) = plot(xform(%s%s(%s)), ''%s%s'');', y_val, suffix, indeces, curr_color, '-');
                end
                try
                        eval(plot_command);
                        if isvar('lstring')
                                legend_str = [legend_str; [' ' lstring{ii}]];
                        else
                                legend_str = [legend_str; sprintf(' %s,%d', suffix(2:end), jj)];
                        end
                catch
                        printf('problem with command %s', plot_command);
                end
		hold on; 
        end
        if eval(sprintf('size(%s, 1) == niters + 1', curr_name))
		orn_index = 'orn_ndx,jj';
                orn_indeces = 'start_ndx:orn_ndx,jj';
	else
		orn_index = 'jj,orn_ndx';
                orn_indeces = 'jj,start_ndx:orn_ndx';
        end
%         if start_ndx == 1
%                 tcomp = sprintf(' - %s(1)', curr_name);
%         end
        % add ornaments
        for jj = 1:Nd
		curr_marker = markers(mod(ii, length(markers)) + 1);
                if plot_ndx == 1
                        plot_command = sprintf('h(ii,2) = plot(sum(%s(%s))%s, xform(%s%s(%s)), ''%s%s'');', curr_name, orn_indeces, tcomp, y_val, suffix, orn_index, curr_color, curr_marker);
                else
                        plot_command = sprintf('h(ii,2) = plot(orn_ndx, xform(%s%s(%s)), ''%s%s'');', y_val, suffix, orn_index, curr_color, curr_marker);
                end
                try
                        eval(plot_command);
                catch
                        printf('problem with command %s', plot_command);
                end
		hold on; 
	end
end
[lh, oh] = legend(h(:,1), legend_str, 'FontSize', FontSize);
for ii = 1:length(exps)
        oh_ndx = length(exps) + (ii - 1)*2 + 1;
        curr_marker = markers(mod(ii, length(markers)) + 1);
        set(oh(oh_ndx), 'Marker', curr_marker);
end
% c = get(l,'Children');
% set(c(1),'Marker','o')
xlabel(xlabels{plot_ndx}, 'FontSize', FontSize)
title(sprintf('%s to x^* over time for [%d %d]', upper(y_val), Nx, Ny))
ylabel('NRMSD to x^* (dB)', 'FontSize', FontSize);
if isvar('plot_axes')
        axis(plot_axes{plot_ndx})
else
        axis tight
end
set(gca,'FontSize', FontSize)
end
