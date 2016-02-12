function [xALP2, CALP2, TALP2, EALP2, ERRALP2] = Ramani_ALP2_wrapper(y, F, S, R, lambda, xinit, niters)
%function [xALP2, CALP2, TALP2, EALP2, ERRALP2] = Ramani_ALP2_wrapper(y, F, S, R, lambda, xinit, niters)
% using the AL-P2 algorithm from Mar. 2011 T-MI paper.
% Sathish Ramani, circa 2011
% 2013-07-14, modified by Jeff Fessler

img = zeros(size(xinit)); % don't know true image
rs = F.arg.Nx;
cs = F.arg.Ny;
Npix = rs*cs;
params.rs = rs;
params.cs = cs;
params.mn = 0; % unknown true min and max vals for experiment
params.mx = 0;
params.Npix = Npix;


smap = S.arg.smaps; 
ncoils = size(smap, 3);
params.ncoils = ncoils;


SP = F.arg.samp; 
lSP = length(find(SP>0)); % # of sampling points

SP3 = repmat(SP, [1 1 ncoils]); 

Data = embed(y, SP3); % Sathish likes to use full matrix version with zeros

	%% Simple iFFT / sum-of-squares recon of raw k-space data
	if ~isvar('recon_sos')
		recon_sos = ifft2(Data);
		recon_sos = sqrt(sum(abs(recon_sos).^2, 3));
		params.recon_sos = recon_sos;
		im(4, recon_sos, 'SoS recon')
	end

if 0
	%% Initial Estimate = SoS + slight perturbation
	if ~isvar('xini')
		randn('state', seedr);
		xini = recon_sos + (randn(rs,cs) + sqrt(-1) * randn(rs,cs))*0.001;
		errini = xini - img;
		errini = sqrt(sum(abs(errini(:)).^2)/Npix)/mx;
		xlabelf('nrmse %g dB', 20*log10(nrms(xini, img)))
		im(6, abs(xini-img), '|SoS error|')
	prompt
	end
end
xini = xinit;

%% Noise decorrelation in data and redefined
% (square root of inverse of noise correlation matrix weighted) sensitivity maps
if ~isvar('smap2')
	recon_F = Npix * ifft2(SP3 .* Data);
	recon_SF = sum(conj(smap) .* recon_F, 3);

	params.smap = smap;
	smap2 = sum(abs(smap).^2, 3); % S'S
%	im(smap2)
end


	params.Prior.PriorType{1} = 'TV'; % Type of penalty for wavelet coef.
	params.Prior.PriorType{2} = 'TV'; % Type of penalty for wavelet coef.

	% Type of sparsifying transform
	params.Operator = 'FD'; % finite differences
        % params.Operator = 'W'; % Redundant wavelet transform
        % params.Operator = 'WFD'; % Redundant wavelet transform & finite diff.

        % Wavelet Options
        params.Wavelet.redundancy = 'undecimated'; % Undecimated wavelet transform using wavelet filters corresponding to standard orthonormal wavelets
        params.Wavelet.wname = 'haar'; % undecimated Haar wavelets
        params.Wavelet.nlev = 2; % wavelet transform with 2 levels
        params.Wavelet.includeApprox = false; % exclude approximation level in the regularizer

%       todo: make it ok in octave
%       Uw = Gwave2('mask', true(rs,cs), 'nlevel', params.Wavelet.nlev)'; % trick

        dwtmode('per', 'nodisp'); % Period boundaries for wavelet implementation
        % params.Wavelet.redundancy = 'none'; % decimatSENSERecon_MRI029_Spiral6_23_2.1_SNR20_OpFD_PriorTV_Lam13.5_variablesed Orthonomal wavelet transform

        % Wavelet filters
        [lod, hid, lor, hir] = wfilters(params.Wavelet.wname);

        % Normalize filters so as to avoid a product with 0.5 during inverse undecimated wavelet transform
        params.Wavelet.lod = lod/sqrt(2);
        params.Wavelet.hid = hid/sqrt(2);
        params.Wavelet.lor = lor/sqrt(2);
        params.Wavelet.hir = hir/sqrt(2);

	% Regularization parameters for the wavelet and TV prior
	params.lambda = lambda;
	params.precon = 0;
	params.betaD = 1e-10;
	params.xinf = zeros(rs,cs);
	params.xinfnorm = 1;
	params.dcosttol = 1e18;
	params.dxtol = 0;

	% Parameters for size of vectors
	[lzALP2 sALP2 eALP2 sr er] = get_Size_AuxVar(params);

	params.AL.lzALP2 = lzALP2;
	params.AL.sALP2 = sALP2;
	params.AL.eALP2 = eALP2;

	params.AL.sr = sr;
	params.AL.er = er;

	RR = compute_RR(params); % FFT of R^T * R
	mxRR = max(RR(:));
	mnRR = min(RR(:));
	params.AL.RR = RR;

	% Other setting
	params.figno = 1;
	params.subplot_img = 5;
	params.subplot_err = 6;
	params.maxitr = niters;
	params.maxitr_in = 1;
	params.dispfig = 1;
	params.dispitr = 1;
	params.dispitrnum = 10;

%% prepare to run AL-P2
if ~isvar('zALP2')
	mxsmap2 = max(smap2(:));
	mnsmap2 = min(smap2(:));
	condSS =  mxsmap2 / mnsmap2; % Condition number of S'S

	kapx = 0.9 * condSS;
	kapu0 = 24; % condition number of (FF + mu*I)
	kapu2 = 12; % condition number of (RR + nu2/n1*I)

	nu2 = (mxsmap2 - mnsmap2 * kapx) / (kapx - 1);
	mu = Npix / (kapu0 - 1); % min eig val of F'F = 0, while max eig val of F'F = N where the DFT matrix F is assumed to be un-normalized 
	nu1 = (kapu2 - 1) * nu2 / (mxRR - mnRR * kapu2);
	nuratio = nu2 / nu1;

	% warn if any of the mu's is negative or less than eps
	if ((mu <= 0) || (nu1 <= 0) || (nu2 <= 0))
		error('At least one of the mu, nu''s is not positive');
	end

	% Inverse of some matrices required for solving sub-problems
	iPpmu = 1 ./ (Npix * SP3 + mu); % Freq. Response of (F'F + mu * I)^-1
	iSpnu2 = 1./(smap2 + nu2); % Inverse of (S'S + nu2 * I)
	iRpnu2nu1 = 1./(RR + nuratio); % Freq. Response of (R'R + nu2/nu1 * I)^-1

	params.AL.mu = mu;
	params.AL.nu1 = nu1;
	params.AL.nu2 = nu2;
	params.AL.iPpmu = iPpmu;
	params.AL.iSpnu2= iSpnu2;
	params.AL.iRpnu2nu1 = iRpnu2nu1;

	dALP2 = zeros(rs, cs, params.AL.lzALP2);
	zALP2 = dALP2;
end

	[xALP2 CALP2 TALP2 EALP2 ERRALP2] = SENSERecon_ALP2(img, ...
		SP3, Data, recon_F, dALP2, zALP2, xini, params);
