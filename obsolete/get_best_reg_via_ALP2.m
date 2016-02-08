%setup;
%addpath(genpath('mai'));
%clear; close all;
%addpath(genpath('..'));

dbug = 0; %simulate baby data
if(dbug)
	nx = 10;
	ny = 12;
else
	nx = 240;
	ny = 200;
% 	nx = 10;
% 	ny = 12;
end
dims = [nx ny];
nc = 4;
if(dbug)
	img = zeros(nx,ny);
	img(nx-8:nx-2,ny-10:ny-3) = 1;
	img(nx-7:nx-5,ny-8:ny-6) = 2;
	sense_maps = zeros(nx,ny,nc);
	sense_maps(1:ceil(nx/2),:,1) = 1;
	sense_maps(floor(nx/2):nx,:,2) = 1;
	sense_maps(:,1:ceil(ny/2),3) = 1;
	sense_maps(:,floor(ny/2):ny,4) = 1;

	sense_maps = zeros(nx,ny,nc);
	sense_maps(1:ceil(2*nx/3),:,1) = 1-rand(size(sense_maps(1:ceil(2*nx/3),:,1)))*0.1 + i*(1-rand(size(sense_maps(1:ceil(2*nx/3),:,1)))*0.1);
	sense_maps(floor(nx/3):nx,:,2) = 1-rand(size(sense_maps(floor(nx/3):nx,:,2)))*0.1 + i*(1-rand(size(sense_maps(floor(nx/3):nx,:,2)))*0.1);
	sense_maps(:,1:ceil(2*ny/3),3) = 1-rand(size(sense_maps(:,1:ceil(2*ny/3),3)))*0.1 + i*(1-rand(size(sense_maps(:,1:ceil(2*ny/3),3)))*0.1);
	%sense_maps(:,floor(ny/3):ny,4) = 1-rand(size(sense_maps(:,floor(ny/3):ny,4)))*0.1 + i*(1-rand(size(sense_maps(:,floor(ny/3):ny,4)))*0.1);
	sense_maps(:,:,4) = 1;
	% dummy smaps
	%sense_maps = ones(nx,ny,nc);

	samp = true(nx,ny);
	samp(5,6) = false;
	samp(3,8) = false;
	samp(8,4) = false;
else
	img = make_sim_image(nx,ny);
%     	img = zeros(nx,ny);
% 	img(nx-8:nx-2,ny-10:ny-3) = 1;
% 	img(nx-7:nx-5,ny-8:ny-6) = 2;
	sense_maps = mri_sensemap_sim('nx',nx,'ny',ny,'ncoil',nc, 'rcoil', 700 );
	% how can I suppress figures from mri_sensemap_sim?
	%figure; im(sense_maps);
	samp = gen_poisson_sampling_pattern('2D',[nx ny], [floor(nx/6) floor(ny/6)],4); % CAUTION: concatenating real poisson disk samplings
	%samp = gen_poisson_sampling_pattern('2D',[nx ny], [floor(nx/3) floor(ny/3)],6); % CAUTION: concatenating real poisson disk samplings
	%keyboard;
end
imgs = repmat(img,[1 1 nc]).*sense_maps;
S = construct_sensefat(sense_maps);
if(dbug)
	imgs_test = reshape(S*img(:),nx,ny,nc);
end




for coil_ndx = 1:nc
	y_full(:,:,coil_ndx) = fft2(imgs(:,:,coil_ndx)); 
	if(dbug)
	    im_test(:,:,coil_ndx) = ifft2(y_full(:,:,coil_ndx));
	    y_full_test(:,:,coil_ndx) = fft2(imgs_test(:,:,coil_ndx));
	end
end
if(dbug)
	y = get_elts(y_full,repmat(samp,[1 1 nc]))';
	ytest3 = get_elts(y_full_test,repmat(samp,[1 1 nc]))';
end
%% construct fatrices and choose parameters
% construct fatrices

% subsampling Fourier encoding matrix F
F = construct_fourierfat(samp,nc);

% sensitivity encoding matrix S
% for coil_ndx = 1:nc
% test_sense_maps(:,:,coil_ndx) = sense_maps(:,:,coil_ndx).';
% end
%S = construct_sensefat(sense_maps);

if(dbug)
	% compare to SoS solution
	zero_fill = reshape(F'*y,nx,ny,nc)/(nx*ny);
	SoS = sqrt(sum(zero_fill.^2,3));
	figure; subplot(2,2,1); im(SoS);

	title('sum of squares solution');

	ytest = F*(S*img(:));
	zero_fill_test = reshape(F'*ytest,nx,ny,nc)/(nx*ny);
	SoS_test = sqrt(sum(zero_fill_test.^2,3));
	subplot(2,2,2); im(SoS_test);
	title('use F and S');

	ytest2 = F*(im_test(:));
	zero_fill_test2 = reshape(F'*ytest2,nx,ny,nc)/(nx*ny);
	SoS_test2 = sqrt(sum(zero_fill_test2.^2,3));
	subplot(2,2,3); im(SoS_test2);
	title('only use F');

	zero_fill_test3 = reshape(F'*ytest3,nx,ny,nc)/(nx*ny);
	SoS_test3 = sqrt(sum(zero_fill_test3.^2,3));
	subplot(2,2,4); im(SoS_test3);
	title('only use S');

	keyboard;
else

	y = F*(S*img(:));
	snr = 40; % specify desired SNR (of sampled values!) in dB
	snr = 80;
	%sig = 10^(-snr/20) * norm(ytrue(samp)) / sqrt(sum(samp(:)));
	sig = 10^(-snr/20) * norm(y) / sqrt(length(y));
	noisey = y + sig*randn(size(y)) + i*sig*randn(size(y));
end
% finite-differencing matrices C1 and C2
[C1, C2] = construct_finite_diff(dims);
%close; % closes window that pops up from using C2sparse
%keyboard;
%C1C1 = full(C1'*C1);
%C2C2 = full(C2'*C2);

% choose parameters


% beta, spatial regularization parameter
%beta = 0.6;
%betas = 0.01:0.01:5;
%betas = 0:0.25:10;
betas=  logspace(0.2,0.5,10);
%beta = 0.52; % or beta = 5+; % for toy case
%beta = 5;
%betas = 1;
%niters = 2000;
niters = 100;

for beta_ndx = 1:length(betas)
beta = betas(beta_ndx);
%beta = 100;
%beta = 30.5386;

% alpha, tradeoff between S
alph = 0.5;
% mu, convergence parameters
mu(1) = 1;
mu(2) = 1;
mu(3) = 4;%mu2
mu(4) = 4;
mu(5) = 4;

%% 

[xhat_P2,xsaved_P2,err_P2,costOrig_P2] = AL_P2_refurb(noisey,F,S,C1,C2,zeros(nx*ny,1),niters,beta,mu,nx,ny,img);
%load ADMM_ALP2_2000iter_beta42_just_xinfs
for iter = 1:niters
	err(iter) = sqrt(mean((abs(col(xsaved_P2(:,:,iter))-img(:)).^2)));
end	


err_per_iter(:,beta_ndx) = err;
beta_xhat(:,:,:,beta_ndx) = xhat_P2;


end

save('best_reg_via_P2.mat','betas','alph','err_per_iter','beta_xhat','img','noisey','F','S','C1','C2','mu','snr');

