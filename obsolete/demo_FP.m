% demo SENSE_fp

dbug = 0; %simulate baby data
small = 0;
if(small)
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
if(small)
	img = zeros(nx,ny);
	img(nx-8:nx-2,ny-10:ny-3) = 1;
	img(nx-7:nx-5,ny-8:ny-6) = 2;
	sense_maps = zeros(nx,ny,nc);
	sense_maps(1:ceil(nx/2),:,1) = 1;
	sense_maps(floor(nx/2):nx,:,2) = 1;
	sense_maps(:,1:ceil(ny/2),3) = 1;
	sense_maps(:,floor(ny/2):ny,4) = 1;

	samp = true(nx,ny);
	samp(5,6) = false;
	samp(3,8) = false;
	samp(8,4) = false;
else
	img = make_sim_image(nx,ny);
%     	img = zeros(nx,ny);
% 	img(nx-8:nx-2,ny-10:ny-3) = 1;
% 	img(nx-7:nx-5,ny-8:ny-6) = 2;
	sense_maps = mri_sensemap_sim('nx',nx,'ny',ny,'ncoil',nc, 'rcoil', 600 );
	% how can I suppress figures from mri_sensemap_sim?
	%figure; im(sense_maps);
	samp = gen_poisson_sampling_pattern('2D',[nx ny], [floor(nx/6) floor(ny/6)],6); % CAUTION: concatenating real poisson disk samplings
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

% subsampling Fourier encoding matrix F
F = construct_fourierfat(samp,nc);

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
	%sig = 10^(-snr/20) * norm(ytrue(samp)) / sqrt(sum(samp(:)));
	sig = 10^(-snr/20) * norm(y) / sqrt(length(y));
	noisey = y + sig*randn(size(y)) + i*sig*randn(size(y));
end
% finite-differencing matrices C1 and C2
[C1, C2] = construct_finite_diff(dims);


% beta, spatial regularization parameter
betas=  logspace(-3,4,100);
betas = 1;
niters = 80;
alphs = [0 0.5 1];
alphs = 0.5;
for beta_ndx = 1:length(betas)
	beta = betas(beta_ndx);

	% alpha, tradeoff between S
	for alph_ndx = 1:length(alphs)
		alph = alphs(alph_ndx);
		%alph = 0;
		%alph = 1;
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
		%%
		%[xhat_FP, xsaved_FP, cost_FP] = SENSE_FP(noisey,F,S,C1,C2,alph,beta,mu,nx,ny,img,niters,1);

		[xhat, xsaved, cost] = SENSE_FP(noisey,F,S,C1,C2,alph,beta,mu,nx,ny,img,niters,0);

%		%load ADMM_ALP2_2000iter_beta42_just_xinfs
%		for iter = 1:niters
%			%err(iter) = norm((abs(xsaved(:,:,iter)-x2000).^2));
%			err(iter) = sqrt(mean((abs(col(xsaved_FP(:,:,iter))-img(:)).^2)));
%		end	
%
%		err_per_iter(:,beta_ndx,alph_ndx) = err;
%		beta_xhat(:,:,:,beta_ndx,alph_ndx) = xhat_FP;

	end
end

figure; im(reshape(xhat,nx,ny));

% use actual AL_P2
%[xhat_P2,xsaved_P2,err_P2,costOrig_P2] = AL_P2_refurb(noisey,F,S,C1,C2,zeros(nx*ny,1),niters,beta,mu,nx,ny,img);

