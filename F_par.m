function F = F_par(Nx, Ny, Nc)
%function F = F_par(Nx, Ny, Nc)
%
%  parallelized fully sampled FFTs over Nc coils 
% abandoned because slower than GsplineF for 8 coils

arg.Nx = Nx;
arg.Ny = Ny;
arg.Nc = Nc;
arg.Nr = arg.Nx*arg.Ny;

if (arg.Nc > 1)
	F = fatrix2('idim', [arg.Nx arg.Ny arg.Nc], 'arg', ...
		arg, 'odim', [arg.Nx arg.Ny arg.Nc], 'forw', ...
		@GsplineF_forw, 'back', @GsplineF_back);
else
	F = fatrix2('idim', [arg.Nx arg.Ny], 'arg', arg, ...
		'odim', [arg.Nx arg.Ny], 'forw', ...
		@GsplineF_forw, 'back', @GsplineF_back);
end

end

% y = G * x
function S = GsplineF_forw(arg, s)

% % Fourier Transform
% S = zeros(arg.Nx,arg.Ny,arg.Nc);
% 	for coil_ndx = 1:arg.Nc
% 		S(:,:,coil_ndx) = fft2(s(:,:,coil_ndx)); 
% 	end

S = zeros(arg.Nx, arg.Ny, arg.Nc);
parfor coil = 1:arg.Nc
	S(:,:,coil) = fft(fft(s(:,:,coil), [], 1), [], 2);
end

end

% x = G' * y
function s = GsplineF_back(arg, S)

% s = zeros(arg.Nx,arg.Ny,arg.Nt,arg.Nc);
%     for coil_ndx = 1:arg.Nc
% 	s(:,:,coil_ndx) = arg.Nr*ifft2(S(:,:,,coil_ndx));
%     end

s = zeros(arg.Nx, arg.Ny, arg.Nc);
parfor coil = 1:arg.Nc
	s(:,:,coil) = ifft(ifft(S(:,:,coil), [], 1), [], 2) * arg.Nr;
end

end
