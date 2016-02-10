function xMFIS = MFISTA_wrapper(Nx, Ny, R, y, xinit, F, S, lambda, niter)
%function xMFIS = MFISTA_wrapper(Nx, Ny, R, y, xinit, F, S, lambda, niter)
% 
% wrapper for MFISTA for getting x_inf for tridiagonal ADMM

	params.Nx = Nx;
	params.Ny = Ny;
	params.dispitr = 1;
	params.dispitrnum = 10;	
	params.dispfig = 0;
	params.Nmask = Nx*Ny;
	params.clim = [];
	params.scale = 1;
	params.figno = 3;
	params.maxitr = 18;
	params.dxtol = 0;
	params.dcosttol = 1e18;
	params.zoomr = 1:Nx;
	params.zoomc = 1:Ny;
	params.formatstringC = '%0.5E';
	params.formatstringT = '%0.3E';
	params.formatstringO = '%0.2E';

	params.ig.mask = true(Nx, Ny);
	params.img = zeros(Nx, Ny);
	
	params.figno = 2;
	params.doMFISTA = 1;
	params.maxitr = niter; % Max # of outer iterations
	params.R = R;

	data = y; 
	xini = xinit;
	A = F*S;
	AWy = A'*(y);
	params.A = A;
	params.AWy = AWy;
	params.W = ones(size(y));
	params.Operator = 'FD';
	params.PriorType = 'l1';
	params.rw = ones(Nx, Ny);%, 2);
	params.lambda = lambda;
	params.MFISTA.nCG = 5;

	params.N = [Nx Ny];

	if exist('mEAWA_MFISTA.mat')
		load('mEAWA_MFISTA.mat')
	else
		params.eigtol = eps; % Matlab epsilon
		params.eigpowermaxitr = 10000;
		params.eigdispitr = 10;	
		mEAWA = get_MaxEigenValue_1by1(params, 'AWA'); 
		save('mEAWA_MFISTA.mat','mEAWA');
	end
	params.mEAWA = mEAWA;
	params.xinf = zeros(Nx,Ny); % to be populated when x_infinity ( minimizer of the cost ) is available
	params.xinfnorm = 1;
	
	[xMFIS CMFIS TFIS l2DFIS RMSEFIS] = runMFISTA(data, AWy, xini, params);
	if norm(xMFIS - xini) == 0
		display('MFISTA did nothing');
		keyboard;
	end
