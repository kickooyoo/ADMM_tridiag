function [xMFIS, CMFIS, TFIS, l2DFIS, RMSEFIS] = MFISTA_wrapper(Nx, Ny, R, y, xinit, F, S, lambda, niter, curr_folder, varargin)
%function [xMFIS, CMFIS, TFIS, l2DFIS, RMSEFIS] = MFISTA_wrapper(Nx, Ny, R, y, xinit, F, S, lambda, niter, curr_folder, varargin)
% 
% wrapper for MFISTA for getting x_inf for tridiagonal ADMM
	force_new_eig = 0;

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
	params.img = zeros(Nx, Ny); % true image unknown
	
	params.figno = 2;
	params.doMFISTA = 1;
	params.maxitr = niter; 
	params.R = R;

	data = y; 
	xini = xinit;
	A = F*S;
	AWy = A'*(y);
	params.A = A;
	params.AWy = AWy;
	params.W = ones(size(y));
	params.Operator = 'AFD';
	params.PriorType = 'l1';
	params.rw = ones(Nx, Ny);
	params.lambda = lambda;
	params.MFISTA.nCG = 5;

	params.N = [Nx Ny];

	mEAWA_fname = [curr_folder '/mEAWA_MFISTA.mat'];
	if exist(mEAWA_fname) && ~force_new_eig
		load(mEAWA_fname)
	else
		params.eigtol = eps; % Matlab epsilon
		params.eigpowermaxitr = 10000;
		params.eigdispitr = 10;	
                tic
		mEAWA = get_MaxEigenValue_1by1(params, 'AWA'); 
                toc
		keyboard
		save(mEAWA_fname,'mEAWA');
 	end
	params.mEAWA = mEAWA;
	params.xinf = zeros(Nx,Ny); 
	params.xinfnorm = 1;

	params = vararg_pair(params, varargin);
	
	[xMFIS, CMFIS, TFIS, l2DFIS, RMSEFIS] = runMFISTA(data, AWy, xini, params);
	if norm(xMFIS - xini) == 0
		display('problem: MFISTA did nothing, check mEAWA and dgrad sign');
		keyboard;
	end
