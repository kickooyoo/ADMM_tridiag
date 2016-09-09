% test tridiag inpaint
inpainting_setup;
niters = 1000;

save_fname = sprintf('inpainting_mat/%s/timing/inpainting_tuning_%s_iters%d_wavelet%d_SNR%d_reduce%1.2d_%strue.mat', obj, machine(1:3), niters, wavelets, SNR, reduce, true_opt);
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

kapps = [2 4 6];
for aa = 1:length(kapps)
	kapp = kapps(aa);
	[x(:,:,aa), xsaved(:,:,:,aa), err(:,aa), cost(:,aa), time(:,aa)] = AL_tridiag_inpaint(y, D, CHW, CVW, ...
		beta, xinit, xtrue, niters, 'betaw', betaw, 'alphw', alphw, 'kapp', kapp);
end

save(save_fname)

if 0 
	figure; plot(cumsum(time), err)
	hold on; plot(cumsum(time_circ), err_circ, 'r')
	figure; subplot(1,2,1); im(x); subplot(1,2,2); im(x - xtrue);
	figure; subplot(1,2,1); im(x_circ); subplot(1,2,2); im(x_circ - xtrue);
end

