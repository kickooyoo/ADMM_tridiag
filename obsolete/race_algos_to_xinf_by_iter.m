
	nx = 240;
	ny = 200;
dims = [nx ny];
nc = 4;

	img = make_sim_image(nx,ny);
	sense_maps = mri_sensemap_sim('nx',nx,'ny',ny,'ncoil',nc, 'rcoil', 700 );
	samp = gen_poisson_sampling_pattern('2D',[nx ny], [floor(nx/6) floor(ny/6)],6); % CAUTION: concatenating real poisson disk samplings
imgs = repmat(img,[1 1 nc]).*sense_maps;
S = construct_sensefat(sense_maps);




for coil_ndx = 1:nc
	y_full(:,:,coil_ndx) = fft2(imgs(:,:,coil_ndx)); 
end

% subsampling Fourier encoding matrix F
F = construct_fourierfat(samp,nc);


	y = F*(S*img(:));
	snr = 40; % specify desired SNR (of sampled values!) in dB
	snr = 80;
	%sig = 10^(-snr/20) * norm(ytrue(samp)) / sqrt(sum(samp(:)));
	sig = 10^(-snr/20) * norm(y) / sqrt(length(y));
	noisey = y + sig*randn(size(y)) + i*sig*randn(size(y));
% finite-differencing matrices C1 and C2
[C1, C2] = construct_finite_diff(dims);

% choose parameters

% beta, spatial regularization parameter
niters = 1000;

beta = 1.85; 

% alpha, tradeoff between S
alph = 0.5;
% mu, convergence parameters
mu(1) = 1;
mu(2) = 1;
mu(3) = 4;%mu2
mu(4) = 4;
mu(5) = 4;
mu(6) = 1;
mu(7) = 1;
mu(8) = 1;
mu(9) = 1;

load xinf_via_ALP2_xinf_only

xinit = zeros(nx*ny,1);

%% 

% ALP2
[xhat_P2,xsaved_P2,err_P2,costOrig_P2] = AL_P2_refurb(noisey,F,S,C1,C2,xinit,niters,beta,mu(1:3),nx,ny,img);

for iter = 1:niters
	dist_to_xinf_P2(iter) = sqrt(sum(abs(col(xsaved_P2(:,:,iter))-xinf(:)).^2)/(nx*ny));
end

% tridiag solver
[xhat_tri,xsaved_tri,cost_tri] = SENSE_ADMM_ALP2(noisey,F,S,C1,C2,alph,beta,mu(1:5),nx,ny,xinit,niters);

for iter = 1:niters
	dist_to_xinf_tri(iter) = sqrt(sum(abs(col(xsaved_tri(:,:,iter))-xinf(:)).^2)/(nx*ny));
end

% "parallel"
[xhat_FP,xsaved_FP,cost_tri] = SENSE_FP(noisey,F,S,C1,C2,alph,beta,mu,nx,ny,xinit,niters);


for iter = 1:niters
	dist_to_xinf_FP(iter) = sqrt(sum(abs(col(xsaved_FP(:,:,iter))-xinf(:)).^2)/(nx*ny));
end


save('RMSD_per_iter.mat','beta','alph','xhat_P2','xsaved_P2','dist_to_xinf_P2','xhat_tri','xsaved_tri','dist_to_xinf_tri','xhat_FP','xsaved_tri','dist_to_xinf_FP','img','noisey','F','S','C1','C2','mu','snr','xinf');

