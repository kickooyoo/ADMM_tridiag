% test tridiag inpaint
inpainting_setup;
niters = 1000;

save_fname = sprintf('inpainting_mat/inpainting_timing_%s_iters%d_wavelet%d_SNR%d_reduce%1.2d.mat', machine(1:3), niters, wavelets, SNR, reduce);
if exist(save_fname, 'file')
	display('file already exists');
	keyboard;
end
sweep_range = 5;
if wavelets
%	betaw = ??;
%	betaw_circ = ??;
	alphw = 1;
	CHW = [CH; betaw * alphw / beta * W];
	CVW = [CV; betaw * (1-alphw) / beta * W];     
	RW = [CHW; CVW];
	RcircW = [Rcirc; betaw_circ / beta_circ * W];
else
	betaw = 0;
	CHW = CH;
	CVW = CV;
	RW = R;
	RcircW = Rcirc;
end
alphas = [0 0.25 5 0.75 1];
for aa = 1:length(alphas)
	alph = alphas(aa);
	[x(:,:,aa), xsaved(:,:,:,aa), err(:,aa), cost(:,aa), time(:,aa)] = AL_tridiag_inpaint(y, D, CHW, CVW, ...
		beta, xinit, xtrue, niters, 'mu', {mu0, mu1, mu2}, 'betaw', betaw, 'alphw', alphw, 'alph', alph);
end

[x_P2, xsave_P2, err_P2, costOrig_P2, time_P2] = AL_P2_inpainting(y, D, RW, ...
	xinit, niters, beta, xtrue, 'mu',  {mu0, mu1} );

[x_circ, xsave_circ, err_circ, costOrig_circ, time_circ] = AL_P2_inpainting(y, D, RcircW, ...
	xinit, niters, beta, xtrue, 'mu', {mu0, mu1});

save(save_fname)

if 0 
	figure; plot(cumsum(time), err)
	hold on; plot(cumsum(time_circ), err_circ, 'r')
	figure; subplot(1,2,1); im(x); subplot(1,2,2); im(x - xtrue);
	figure; subplot(1,2,1); im(x_circ); subplot(1,2,2); im(x_circ - xtrue);
end

