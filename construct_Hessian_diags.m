function [subCC, subCCT, diagCC, diagCCT] = construct_Hessian_diags(mu0, mu1, mu3, mu4, Nx, Ny, beta, varargin)

arg.betaw = 0;
arg.alphw = 0.5;
arg = vararg_pair(arg, varargin);
if ~isscalar(mu3)
        mu3 = reshape(mu3, Nx, Ny);
end
if ~isscalar(mu4)
        mu4 = reshape(mu4, Nx, Ny);
end

% pass tridiag of C'C into mex
Wconst = arg.betaw * arg.alphw / beta;
WVconst = arg.betaw * (1-arg.alphw) / beta;
if ~isscalar(mu0) || ~isscalar(mu1);
	subCCT = single(-mu1(:,1:end-1).' .* ones(Ny - 1, Nx)); % tranpose dims 
	subCC = single(-mu0(1:end-1,:) .* ones(Nx - 1, Ny));
else
	subCCT = single(-mu1 * ones(Ny - 1, Nx)); % tranpose dims 
	subCC = single(-mu0 * ones(Nx - 1, Ny));
end
if ~isscalar(mu3) || ~isscalar(mu4)
        diagCCT = single(mu1.' .* (cat(1, ones(1, Nx), 2*ones(Ny-2, Nx), ones(1, Nx)) + WVconst.^2) + ...
                mu3.' + mu1.' .* WVconst);
        diagCC = single(mu0 .* (cat(1, ones(1, Ny), 2*ones(Nx-2, Ny), ones(1, Ny)) + Wconst.^2) + ...
                mu4 + mu0 .* Wconst);
else
        diagCCT = single(mu1 .* (cat(1, ones(1, Nx), 2*ones(Ny-2, Nx), ones(1, Nx)) + WVconst.^2) + ...
                mu3 + mu1 .* WVconst);
        diagCC = single(mu0 .* (cat(1, ones(1, Ny), 2*ones(Nx-2, Ny), ones(1, Ny)) + Wconst.^2) + ...
                mu4 + mu0 .* Wconst);
end
if any(imag(cat(1, col(subCC), col(subCCT), col(diagCC), col(diagCCT))))
	display('invalid diagonal values, some are imaginary');
	keyboard;
end