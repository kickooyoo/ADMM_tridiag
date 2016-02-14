% setup tridiag exp
truncate = 0;
wavelets = 0;

if ~isvar('mapped_im')
	fn = [home_path 'Documents/data/2010-07-06-fessler-3d/raw/p23040-3dspgr-8ch.7'];

	nc = 8;
	[im4d,d] = recon3dft(fn,nc);
	[body_coil_images,d] = recon3dft([home_path 'Documents/data/2010-07-06-fessler-3d/raw/p24064-3dspgr-body.7'],1);

	% here we clip the images to see affect of factors of 2 on the FFT
	if truncate
		im4d = im4d(3:end-2,3:end-2,:,:); % size orig [256 144], now [254 140]
	end

	% slice selection
	% sathish used 38, Muckley and Dan complained, so switch to 67
	slice = 67;
	slice = 38;
	mapped_im = squeeze(im4d(:,:,slice,:));
	Nx = size(im4d,1);
	Ny = size(im4d,2);
	dims = [Nx Ny];

	% need to estimate sense maps for other slices
	%load(sprintf('%sDocuments/mai_code/static_SENSE_splitting/saved_results/experimental/smap_est_slice_%d.mat', home_path, slice));
	load(sprintf('%sDocuments/data/2010-07-06-fessler-3d/slice38/ramani/Smaps%d.mat', home_path, slice), 'Smap_QPWLS'); S_est = Smap_QPWLS;
	if truncate
		sense_maps = S_est(3:end-2,3:end-2,:);
	else
		sense_maps = S_est;
	end
	clear S_est;
end

if ~isvar('samp')
	% generate sampling pattern
	if 1
		if exist('PD_sampling_256x144_R6_center8.mat')
			load('PD_sampling_256x144_R6_center8.mat','samp');
		else
			reduction = 6;
			params.Nx = Nx;
			params.Ny = Ny;
			params.Nf = 1;
			params.h = 1;
			params.R = round(Nx*Ny/reduction);
			samp = logical(gen_new_sampling_pattern(params, 'all_kspace', 0, 'ellipse', 0, 'ncenter', 8, 'vardensity', 0));
			save('PD_sampling_256x144_R6_center8.mat','samp');
		end
	else
		load('al-p2/PDiskR256x256L1R2.252.25B0.8R1D1pctg20.2972.mat', 'SP', 'sampname');
		if size(SP, 1) > Nx
			SP = SP(1:Nx,:);
		elseif size(SP, 1) < Nx
			display('need bigger Poisson disk sampling pattern');
			keyboard;
		end
		if size(SP, 2) > Ny
			SP = SP(:,1:Ny);
		elseif size(SP, 2) < Ny
			display('need bigger Poisson disk sampling pattern');
			keyboard;
		end
		Ncentx = 8; % sample a window of 2*N around DC along x
		Ncenty = 8; % sample a window of 2*N around DC along y

		samp = logical(coverDC_SamplingMask(SP, Ncentx, Ncenty));

	end
end

% construct fatrices
S = GsplineS(sense_maps, 1);
F = GsplineF(Nx, Ny, 1, nc, 'samp', samp);
[CH, CV] = construct_finite_diff(dims); 
R = [CH; CV];

% generate data
y = F * mapped_im(:);
SNR = 40;
sig = 10^(-SNR/20) * norm(y) / sqrt(length(y));
rng(0, 'twister')
y_noise = y + sig*randn(size(y)) + 1i*sig*randn(size(y));


% initialize with SoS zero-fill solution
center_samp = logical(coverDC_SamplingMask(zeros(Nx, Ny), 16, 16));
F_center = GsplineF(Nx, Ny, 1, nc, 'samp', center_samp);
y_center = F_center * mapped_im(:);
y_center_noise = y_center + sig*randn(size(y_center)) + 1i*sig*randn(size(y_center));
zero_fill = reshape(F_center'*y_center_noise, Nx, Ny, nc)/(Nx * Ny);
SoS = sqrt(sum(abs(zero_fill).^2,3));
% [xinit, scale] = ir_wls_init_scale(F*S, y_noise, SoS);
%SoS_compensate = sum(conj(sense_maps).*mapped_im,3)./(sum(abs(sense_maps).^2,3));
xinit = SoS;
if 0
        [xinit, scale] = ir_wls_init_scale(A, y_center_noise, SoS);
        [xinit_MFIS, scale_MFIS] = ir_wls_init_scale(A, y_center_noise, xMFIS);
        figure; im(cat(1, cat(2, SoS, xinit), cat(2, xMFIS, xinit_MFIS)))
        title(sprintf('SoS scale: %.2d + %.2di, MFISTA scale: %.2d + %.2di',real(scale),imag(scale),real(scale_MFIS),imag(scale_MFIS)))
        [~, testscale] = ir_wls_init_scale(1, xMFIS, SoS)
        % testscale = 0.1600 - 0.4581i
end

if slice == 67
	mask = generate_mask('slice67',1,Nx,Ny);
	if truncate
		mask = mask(3:end-2,3:end-2);
	end
elseif slice == 38
	mask = generate_mask('slice38',1,Nx,Ny);
else
	mask = true(Nx, Ny);
end

% parameters
if slice == 67
	beta = 2^19;
elseif slice == 38
	beta = 2^24; % for slice 38 and l2b = 12 samp
	beta = 2^28; % for slice 38 and l2b = 16 samp
	beta = 2^25; % for slice 38 and sathish samp and sathish samp
	beta = 2^20; % for slice 38 and new samp R=6, sathish smap
else
	keyboard;
end
% tradeoff between S and x, between [0 1]
alph = 0.5;
% convergence parameters
Rcirc = Cdiffs([Nx Ny],'offsets', [1 Nx], 'type_diff','circshift');
if wavelets
        W = Godwt1(true(Nx, Ny));
        betaw = beta;
        alphw = 0.5;
        RcircW = [Rcirc; W];
        CHW = [CH; betaw * alphw / beta * W];
        CVW = [CV; betaw * (1-alphw) / beta * W];     
	RW = [CHW; CVW];
end
