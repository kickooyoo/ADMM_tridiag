% test tridiag inpaint
inpainting_setup;
niters = 2000;

save_fname = sprintf('inpainting_mat/%s/timing/inpainting_tuning_%s_iters%d_wavelet%d_SNR%d_reduce%1.2d_%strue.mat', obj, machine(1:3), niters, wavelets, SNR, reduce, true_opt);
if exist(save_fname, 'file')
	display('file already exists');
	keyboard;
end

kapps = [2 3 4 6 8 12];
poss = [0.01 0.1 1];
for aa = 1:length(kapps)
	for pp = 1:length(poss)
		pos = poss(pp);
		kapp = kapps(aa);
		[x(:,:,aa,pp), ~, err(:,aa,pp), cost(:,aa,pp), time(:,aa,pp)] = AL_tridiag_inpaint(y, D, CHW, CVW, ...
			beta, xinit, xtrue, niters, 'betaw', betaw, 'alphw', alphw, 'kapp', kapp, 'pos', pos);
	end
end

save(save_fname)

if 0 
	figure; plot(cumsum(time), err)
	hold on; plot(cumsum(time_circ), err_circ, 'r')
	figure; subplot(1,2,1); im(x); subplot(1,2,2); im(x - xtrue);
	figure; subplot(1,2,1); im(x_circ); subplot(1,2,2); im(x_circ - xtrue);
end

