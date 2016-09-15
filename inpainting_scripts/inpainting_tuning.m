% test tridiag inpaint
niters = 2000;

save_fname = sprintf('inpainting_mat/%s/timing/inpainting_tuning_%s_iters%d_wavelet%d_SNR%d_reduce%1.2d_%strue.mat', obj, machine(1:3), niters, wavelets, SNR, reduce, true_opt);
if exist(save_fname, 'file')
	display('file already exists');
	keyboard;
end

kapps = [6 7 8 10 12];
poss = [0.005 0.01 0.05];
for aa = 1:length(kapps)
	for pp = 1:length(poss)
		pos = poss(pp);
		kapp = kapps(aa);
		[x(:,:,aa,pp), ~, err(:,aa,pp), cost(:,aa,pp), time(:,aa,pp)] = AL_tridiag_inpaint(y, D, CHW, CVW, ...
			beta, xinit, xtrue, niters, 'betaw', betaw, 'alphw', alphw, 'kapp', kapp, 'pos', pos);
	end
	save(save_fname)
end


if 0 
	figure; plot(cumsum(time), err)
	hold on; plot(cumsum(time_circ), err_circ, 'r')
	figure; subplot(1,2,1); im(x); subplot(1,2,2); im(x - xtrue);
	figure; subplot(1,2,1); im(x_circ); subplot(1,2,2); im(x_circ - xtrue);
end

