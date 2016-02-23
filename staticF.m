function F = staticF(Nx, Ny, Nc, varargin)
%function F = GsplineF(Nx, Ny,, Nc, varargin)
%
%if nargin < 1, help(mfilename), error(mfilename), end
%if nargin == 1 && streq(varargin{1}, 'test'), Gdft_test, return, end
%
% object f should be Nx x Ny x Nt
% varargin: samp for subsampling (dynamic)

arg.Nx = Nx;
arg.Ny = Ny;
arg.Nc = Nc;

arg.Nr = arg.Nx*arg.Ny;
arg.Ns = arg.Nr;
arg.samp = []; 
arg = vararg_pair(arg, varargin);

if ~isempty(arg.samp) && ( (size(arg.samp,1) ~= Nx) || (size(arg.samp,2) ~= Ny) )
        display(sprintf('samp size [%d %d] does not match given dims:[%d %d]', size(arg.samp,1), size(arg.samp,2), Nx, Ny));
        keyboard
end
if ~isempty(arg.samp) && ~islogical(arg.samp)
        display('samp must be logical');
        keyboard;
end
if ~isempty(arg.samp)
        arg.Ns = sum(col(arg.samp));
end

if (arg.Nc > 1)
        if isempty(arg.samp)
                F = fatrix2('idim', [arg.Nx arg.Ny arg.Nc], 'arg', ...
                        arg, 'odim', [arg.Nx arg.Ny arg.Nc], 'forw', ...
                        @staticF_forw, 'back', @staticF_back);
        else
                F = fatrix2('idim', [arg.Nx arg.Ny arg.Nc], 'arg', ...
                        arg, 'odim', [arg.Ns arg.Nc], 'forw', ...
                        @staticF_forw, 'back', @staticF_back);
                
        end
else
        if isempty(arg.samp)
                F = fatrix2('idim', [arg.Nx arg.Ny], 'arg', arg, ...
                        'odim', [arg.Nx arg.Ny], 'forw', ...
                        @staticF_forw, 'back', @staticF_back);
        else
                F = fatrix2('idim', [arg.Nx arg.Ny], 'arg', arg, ...
                        'odim', [arg.Ns], 'forw', @staticF_forw, ...
                        'back', @staticF_back);
                
        end
end

end

% y = G * x
function S = staticF_forw(arg,s)

S = fft(fft(s,[],1),[],2);
if ~isempty(arg.samp)
        S = S(repmat(arg.samp, [1 1 arg.Nc]));
end

end

% x = G' * y
function s = staticF_back(arg,S)

if ~isempty(arg.samp)
        S = embed(col(S), repmat(arg.samp, [1 1 arg.Nc]));
end
s = ifft(ifft(S,[],1),[],2)*arg.Nr;

end
