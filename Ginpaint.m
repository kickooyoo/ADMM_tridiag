function G = Ginpaint(samp, varargin)
% function G = Ginpaint(samp, varargin)
% 

arg.Nx = size(samp, 1);
arg.Ny = size(samp, 2);
arg.samp = samp;
arg.Nsamp = sum(col(samp));
arg.oldfat = 0; % use fatrix, not fatrix2
arg.mask = true(size(samp)); % to induce imask
arg = vararg_pair(arg, varargin);

if arg.oldfat
%	G = Gmatrix(diag(samp(:)), 'idim', [arg.Nx arg.Ny], 'odim', arg.Nsamp); % too much memory!
	G = Fatrix('idim', [arg.Nx arg.Ny] ,'arg', arg,'odim', ...
		arg.Nsamp, 'forw', @Ginpaint_forw, 'back', @Ginpaint_back);
else
	G = fatrix2('idim', [arg.Nx arg.Ny] ,'arg', arg,'odim', ...
		arg.Nsamp, 'forw', @Ginpaint_forw, 'back', @Ginpaint_back);
end
end

% y = G * x
function y = Ginpaint_forw(arg,x)

y = masker(x, arg.samp);

end

% x = G' * y
function x = Ginpaint_back(arg,y)

x = embed(y, arg.samp);
end
