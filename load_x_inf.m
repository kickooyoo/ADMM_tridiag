function xinf = load_x_inf(slice, beta, curr_folder, slice_str, varargin)
%function xinf = load_x_inf(slice, beta, curr_folder, slice_str, varargin)
%
% varargin: 
%	truncate (boolean)
%	method (string) 'MFISTA' (default) or 'tridiag'

arg.truncate = 0;
arg.method = 'MFISTA';
arg = vararg_pair(arg, varargin);

if strcmp(arg.method, 'tridiag')
	% load xinfs
	fname = sprintf('%s/x_tri_inf_%d_beta%.*d.mat', curr_folder, slice_str, 3, beta);
	if exist(fname, 'file')
		load(fname, 'x_tri_inf');
		if arg.truncate
			x_tri_inf = reshape(x_tri_inf, 256, 144);
			xinf = x_tri_inf(3:end-2, 3:end-2);
		else
			xinf = x_tri_inf;
		end
		display('using tridiag ADMM as true!')
	else
		display(sprintf('unable to load xinf from file %s', fname));
		xinf = 0;
	end
elseif strcmp(arg.method, 'MFISTA')
	fname = sprintf('%s/x_MFISTA_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta);
	if exist(fname, 'file')
		load(fname, 'x*MFIS*');
		if isvar('xMFIS')
			x_MFISTA = xMFIS;
		end
		if arg.truncate
			xMFIS = reshape(x_MFISTA, 256, 144);
			xinf = xMFIS(3:end-2, 3:end-2);
		else
			xinf = x_MFISTA;
		end
	else
		display(sprintf('unable to load xinf from file %s', fname));
		xinf = 0;
	end
else
	xinf = 0;
	display('unknown method of x inf');
	keyboard;
end

