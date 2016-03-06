function xinf = load_x_inf(slice, beta, curr_folder, slice_str, varargin)
%function xinf = load_x_inf(slice, beta, curr_folder, slice_str, varargin)
%
% varargin: 
%	truncate (boolean)
%	method (string) 'avg' (default), 'tridiag', or 'MFISTA'

arg.truncate = 0;
arg.method = 'avg';%'MFISTA';
arg = vararg_pair(arg, varargin);

if strcmp(arg.method, 'tridiag') || strcmp(arg.method, 'avg')
	% load xinfs
	fname = sprintf('%s/x_tri_inf_%d_beta%.*d.mat', curr_folder, slice_str, 3, beta);
	if exist(fname, 'file')
		load(fname, 'x_tri_inf');
		if arg.truncate
			x_tri_inf = reshape(x_tri_inf, 256, 144);
			x_tri_inf = x_tri_inf(3:end-2, 3:end-2);
		end
		xinf = x_tri_inf;
		display('using tridiag ADMM as true!')
	else
		display(sprintf('unable to load xinf from file %s', fname));
		xinf = 0;
	end
end
if strcmp(arg.method, 'MFISTA') || strcmp(arg.method, 'avg')
	fname = sprintf('%s/x_MFISTA_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta);
	if exist(fname, 'file')
		load(fname, 'x*MFIS*');
		if isvar('xMFIS')
			x_MFISTA_inf = xMFIS;
		endi
		if arg.truncate
			x_MFISTA_inf = reshape(x_MFISTA, 256, 144);
			x_MFISTA_inf = x_MFISTA_inf(3:end-2, 3:end-2);
		end
		xinf = x_MFISTA_inf;
	else
		display(sprintf('unable to load xinf from file %s', fname));
		xinf = 0;
	end
end
if strcmp(arg.method, 'avg')
	if ~isvar('x_tri_inf') || ~isvar('x_MFISTA_inf')
		display('not able to load one of xinfs for average');
		keyboard;
	end
	xinf = mean(cat(3, x_tri_inf, x_MFISTA_inf), 3);
end
if ~isvar('xinf')
	xinf = 0;
	display('unknown method of x inf');
	keyboard;
end

