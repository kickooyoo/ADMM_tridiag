function S = staticS(sense_maps)
%function S = GsplineS(sense_maps)
% 

assert(ndims_ns(sense_maps) == 3, 'invalid dimensions for sense_maps');

arg.Nx = size(sense_maps,1);
arg.Ny = size(sense_maps,2);
arg.Nr = arg.Nx*arg.Ny;
arg.Nc = size(sense_maps,3);
arg.mask = true(arg.Nx, arg.Ny);
arg.smaps = sense_maps;

S = fatrix2('idim', [arg.Nx arg.Ny] ,'arg',arg,'odim', ...
        [arg.Nx arg.Ny arg.Nc], 'forw', @staticS_forw, ...
        'back', @staticS_back, 'imask', arg.mask);
end

% y = G * x
function s = staticS_forw(arg,x)

x_rep = repmat(x,[1 1 arg.Nc]);
s = x_rep.*arg.smaps;

end

% x = G' * y
function x = staticS_back(arg,s)

prod_rep = conj(arg.smaps).*s;
x = sum(prod_rep,3);

end
