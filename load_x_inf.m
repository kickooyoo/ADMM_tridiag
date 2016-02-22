function x_inf = load_xinf(slice, beta, varargin)
%function x_inf = load_xinf(slice, beta, varargin)
%
% varargin: 
%	truncate (boolean)
%	method (string) 'MFISTA' (default) or 'tridiag'

arg.truncate = 0;
arg.method = 'MFISTA';
arg = vararg_pair(arg, varargin);

if slice == 0
	slice_str = 'sim';
else
	slice_str = sprintf('slice%d', slice);
end

if strcmp(arg.method, 'tridiag')
	% load xinfs
	load(sprintf('./reviv/curr/x_tri_inf_%d_beta%.*d.mat', slice_str, 3, beta), 'x_tri_inf');
	if arg.truncate
		x_tri_inf = reshape(x_tri_inf, 256, 144);
		x_inf = x_tri_inf(3:end-2, 3:end-2);
	else
		x_inf = x_tri_inf;
	end
	display('using tridiag ADMM as true!')
elseif strcmp(arg.method, 'MFISTA')
	load(sprintf('./reviv/curr/x_MFISTA_inf_%s_beta%.*d.mat', slice_str, 3, beta), 'x*MFIS*');
	if isvar('xMFIS')
		x_MFISTA = xMFIS;
	end
	if arg.truncate
		xMFIS = reshape(x_MFISTA, 256, 144);
		x_inf = xMFIS(3:end-2, 3:end-2);
	else
		x_inf = x_MFISTA;
	end
else
	x_inf = 0;
	display('unknown method of x inf');
	keyboard;
end

