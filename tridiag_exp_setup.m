% setup tridiag exp
truncate = 0;
wavelets = 0;

if ~isvar('mapped_im')
	fn = [home 'Documents/data/2010-07-06-fessler-3d/raw/p23040-3dspgr-8ch.7'];

	nc = 8;
	[im4d,d] = recon3dft(fn,nc);
	[body_coil_images,d] = recon3dft([home 'Documents/data/2010-07-06-fessler-3d/raw/p24064-3dspgr-body.7'],1);

	% here we clip the images to see affect of factors of 2 on the FFT
	if truncate
		im4d = im4d(3:end-2,3:end-2,:,:); % size orig [256 144], now [254 140]
	end

	% slice selection
	% sathish used 38, Muckley and Dan complained, so switch to 67
	slice = 67;
	mapped_im = squeeze(im4d(:,:,slice,:));
	Nx = size(im4d,1);
	Ny = size(im4d,2);
	dims = [Nx Ny];

	% need to estimate sense maps for other slices
	% for slice 67
	load([home 'Documents/mai_code/static_SENSE_splitting/saved_results/experimental/smap_est_slice_67.mat']);
	if truncate
		sense_maps = S_est(3:end-2,3:end-2,:);
	else
		sense_maps = S_est;
	end
	clear S_est;
end

if ~isvar('samp')
	% generate sampling pattern
	reduction = 6;
	params.Nx = Nx;
	params.Ny = Ny;
	params.Nf = 1;
	params.h = 1;
	params.R = round(Nx*Ny/reduction);
	samp = logical(gen_new_sampling_pattern(params, 'all_kspace', 0));
end

% construct fatrices
S = GsplineS(sense_maps, 1);
F = GsplineF(Nx, Ny, 1, nc, 'samp', samp);
[CH, CV] = construct_finite_diff(dims); 
R = [CH; CV];

% generate data
y = F*(mapped_im(:));
SNR = 40;
sig = 10^(-SNR/20) * norm(y) / sqrt(length(y));
rng(0)
y_noise = y + sig*randn(size(y)) + 1i*sig*randn(size(y));


% initialize with SoS zero-fill solution
zero_fill = reshape(F'*y,Nx,Ny,nc)/(Nx*Ny);
SoS = sqrt(sum(abs(zero_fill).^2,3));
[xinit, scale] = ir_wls_init_scale(F*S, y_noise, SoS);
%SoS_compensate = sum(conj(sense_maps).*mapped_im,3)./(sum(abs(sense_maps).^2,3));

mask = generate_mask('slice67',1,Nx,Ny);
if truncate
	mask = mask(3:end-2,3:end-2);
end

% parameters
beta = 2^19;
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
