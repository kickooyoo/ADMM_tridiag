function cost = MFISTA_cost_wrapper(Nx, Ny, R, y, xinit, F, S, lambda)
%function cost = MFISTA_cost_wrapper(Nx, Ny, R, y, xinit, F, S, lambda)
% 
% wrapper for MFISTA's compute_Cost_CT2D to compare against tridiag_ADMM_cost 

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
	
	cost = compute_Cost_CT2D(data, xinit, params);

