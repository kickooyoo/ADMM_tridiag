% setup tridiag exp
truncate = 0;
wavelets = 1;

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

% generate sampling pattern
reduction = 6;
params.Nx = Nx;
params.Ny = Ny;
params.Nf = 1;
params.h = 1;
params.R = round(Nx*Ny/reduction);
samp = logical(gen_new_sampling_pattern(params, 'all_kspace', 0));

% construct fatrices
S = GsplineS(sense_maps, 1);
F = GsplineF(Nx, Ny, 1, nc, 'samp', samp);
[CH, CV] = construct_finite_diff(dims); !
R = [CH; CV];
if wavelets
        W = Godwt1(true(Nx, Ny));
        betaw = beta;
        alphw = 0.5;
        Rcirc = Cdiffs([Nx Ny],'type_diff','circshift');
        RcircW = [Rcirc; W];
end

% generate data
y = F*(mapped_im(:));
% sig = ??
y_noise = y + sig*randn(size(y)) + 1i*sig*randn(size(y));


% initialize with SoS zero-fill solution
zero_fill = reshape(F'*y,Nx,Ny,nc)/(Nx*Ny);
SoS = sqrt(sum(zero_fill.^2,3));
SoS_compensate = sum(conj(sense_maps).*mapped_im,3)./(sum(abs(sense_maps).^2,3));

% mask = generate_mask('slice67',1,Nx,Ny);

% parameters
beta = 2^20.1;
% tradeoff between S and x, between [0 1]
alph = 0.5;
% convergence parameters
mu = ones(1,5);
if wavelets
        CHW = [CH; betaw * alphw / beta * W];
        CVW = [CV; betaw * (1-alphw) / beta * W];     
end
RW = [CHW; CVW];