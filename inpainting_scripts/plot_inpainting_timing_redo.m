db_path = '~/Dropbox/fessler/writing/thesis/figures/';
FontSize = 18;
if ~isvar('colors')
        colors = 'cmgkbry';
end
if ~isvar('markers')
        markers = 'xo+*sd.^v><ph';
end
if ~isvar('ystr')
        ystr = 'NRMSE to x_{true} (dB)';
end
if ~isvar('time_MFIS') && isvar('time_FIS')
        time_MFIS = time_FIS;
        
end

clear *circ
clear time_FIS;
short_slice_str = 'inpaint';
order = [2 3 1 4];
x = x(:,:,1);
err = err(:,1);
time = time(:,1);
% lstring = {'AL-tridiag-inpaint'; 'AL-P2-NC-inpaint'; 'AL-P2-inpaint'};
MarkerSize = 10;
LineWidth = 2;
exps = who('time*');
start_ndx = 1;
end_ndx = 200;%niters+1;
if ~isvar('orn_ndx')
        orn_ndx = round(end_ndx/4);
end
xform = @(x) 20*log10(x);

if ~isvar('FontSize')
	FontSize = 12;
end
if ~isvar('order')
        order = 1:length(exps);
end
y_val = 'err';%'nrmsd'; %'nrmsd'; %'costOrig'
xlabels = {'wall time (s)'; 'iteration number'};
for plot_ndx = 1:2
fhandle = figure('Position', [100 100 400 300]); hold on; 
legend_str = {};
for kk = 1:length(exps)
        ii = order(kk);
	curr_name = exps{ii};
	suffix_ndx = strfind(curr_name, '_');
	suffix = curr_name(suffix_ndx : end);
        if eval(sprintf('size(%s,1) == niters + 1', curr_name))
                eval(sprintf('%s = %s.'';', curr_name, curr_name));
                eval(sprintf('%s%s = %s%s.'';', y_val, suffix, y_val, suffix));
        end
 	if eval(sprintf('prod(size(%s)) ~= niters +1', curr_name))
 		Nd = eval(sprintf('ceil(prod(size(%s))/(niters + 1))', curr_name));
 	else
		Nd = 1;
 	end
	if eval(sprintf('size(%s, 1) == niters + 1', curr_name))
		indices = 'start_ndx:end_ndx,jj';
	else
		indices = 'jj,start_ndx:end_ndx';
        end
        if start_ndx == 1
                tcomp = sprintf(' - %s(1)', curr_name);
        end
	curr_color = colors(mod(kk, length(colors)) + 1);
	for jj = 1:Nd
                if plot_ndx == 1
                        plot_command = sprintf('h(kk,1) = plot(cumsum(%s(%s))%s, xform(%s%s(%s)), ''%s%s'', ''LineWidth'', %d);', curr_name, indices, tcomp, y_val, suffix, indices, curr_color, '-', LineWidth);
                else
                        plot_command = sprintf('h(kk,1) = plot(xform(%s%s(%s)), ''%s%s'', ''LineWidth'', %d);', y_val, suffix, indices, curr_color, '-', LineWidth);
                end
                try
                        eval(plot_command);
                        if isvar('lstring')
                                legend_str = [legend_str; [' ' lstring{ii}]];
                        else
                                legend_str = [legend_str; sprintf(' %s,%d', suffix(2:end), jj)]
                        end
                catch
                        printf('problem with command %s', plot_command);
                end
		hold on; 
        end
        if eval(sprintf('size(%s, 1) == niters + 1', curr_name))
		orn_index = 'orn_ndx,jj';
                orn_indices = 'start_ndx:orn_ndx,jj';
	else
		orn_index = 'jj,orn_ndx';
                orn_indices = 'jj,start_ndx:orn_ndx';
        end
%         if start_ndx == 1
%                 tcomp = sprintf(' - %s(1)', curr_name);
%         end
        % add ornaments
        for jj = 1:Nd
		curr_marker = markers(mod(kk, length(markers)) + 1);
                if plot_ndx == 1
                        plot_command = sprintf('h(kk,2) = plot(sum(%s(%s))%s, xform(%s%s(%s)), ''%s%s'', ''MarkerFaceColor'', ''%s'', ''MarkerSize'', MarkerSize);', curr_name, orn_indices, tcomp, y_val, suffix, orn_index, curr_color, curr_marker, curr_color);
                else
                        plot_command = sprintf('h(kk,2) = plot(orn_ndx, xform(%s%s(%s)), ''%s%s'', ''MarkerFaceColor'', ''%s'', ''MarkerSize'', MarkerSize);', y_val, suffix, orn_index, curr_color, curr_marker, curr_color);
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
% add ornaments to lines
for kk = 1:length(exps)
        ii = order(kk);
        oh_ndx = length(exps) + (ii - 1)*2 + 1;
        curr_marker = markers(mod(ii, length(markers)) + 1);
        curr_color = colors(mod(ii, length(colors)) + 1);
        set(oh(oh_ndx), 'Marker', curr_marker, 'MarkerFaceColor', curr_color, 'MarkerSize', MarkerSize);
end
% c = get(l,'Children');
% set(c(1),'Marker','o')
xlabel(xlabels{plot_ndx}, 'FontSize', FontSize)
title(sprintf('%s to %s over time for [%d %d]', upper(y_val),'x^{(\infty)}',  Nx, Ny))
ylabel(ystr, 'FontSize', FontSize);
if isvar('plot_axes')
        axis(plot_axes{plot_ndx})
else
        axis tight
end
set(gca,'FontSize', FontSize)
end
%%

title('');
save_im(db_path, sprintf('figs/%s_iter_%s', short_slice_str, machine(1:3)))
close;
title('');
xlim([0 50])
save_im(db_path, sprintf('figs/%s_timing_%s', short_slice_str, machine(1:3)))
