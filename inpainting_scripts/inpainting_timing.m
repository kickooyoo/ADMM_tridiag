% test tridiag inpaint
if ~isvar('niters')
	niters = 20000;
end
if ~isvar('str_mod')
	str_mod = '';
end
if do_alph
	alphas = 0:0.1:1;
	save_fname = sprintf('inpainting_mat/%s/timing/inpainting_timing_%s_iters%d_wavelet%d_SNR%d_reduce%1.2d_%strue_retunedmu%s.mat', obj, machine(1:3), niters, wavelets, SNR, reduce, true_opt, str_mod);
else
	alphas = 0.5;
	save_fname = sprintf('inpainting_mat/%s/timing/inpainting_timing_%s_iters%d_wavelet%d_SNR%d_reduce%1.2d_%strue_%1.1dalph_tunedmu%s.mat', obj, machine(1:3), niters, wavelets, SNR, reduce, true_opt, alphas, str_mod);
end
if exist(save_fname, 'file')
	display('file already exists');
	keyboard;
end

aaa = find(alphas == 0.5);
[x(:,:,aaa), ~, err(:,aaa), cost(:,aaa), time(:,aaa)] = AL_tridiag_inpaint(y, D, CHW, CVW, ...
		beta, xinit, xtrue, niters, 'betaw', betaw, 'alphw', alphw, 'alph', alph, 'kapp', 7, 'pos', 0.01);

[x_P2, ~, err_P2, costOrig_P2, time_P2] = AL_P2_inpainting(y, D, RW, ...
	xinit, niters, beta, xtrue, 'alphw', alphw, 'betaw', betaw);

[x_circ, ~, err_circ, costOrig_circ, time_circ] = AL_P2_inpainting(y, D, RcircW, ...
	xinit, niters, beta, xtrue);
[x_MFIS, C_MFIS, time_MFIS, err_MFIS, ~] = MFISTA_inpainting_wrapper(Nx, Ny, RW, y, xinit, D, beta, niters, curr_folder, slice_str, 'xinf', xtrue, 'xinfnorm', norm(col(xtrue),2));
save(save_fname)

[x_ADMM, ~, err_ADMM, cost_ADMM, time_ADMM] = ADMM_tridiag_inpaint(y, D, CHW, CVW, beta, xinit, xtrue, niters, 'betaw', betaw, 'alphw', alphw, 'alph', alph, 'kapp', 7, 'pos', 0.01);
save(save_fname)
send_mai_text('done with inpainting timing except alph')
for aa = setdiff(1:length(alphas), aaa)
	alph = alphas(aa);
	[x(:,:,aa), ~, err(:,aa), cost(:,aa), time(:,aa)] = AL_tridiag_inpaint(y, D, CHW, CVW, ...
		beta, xinit, xtrue, niters, 'betaw', betaw, 'alphw', alphw, 'alph', alph, 'kapp', 7, 'pos', 0.01);
end
save(save_fname)
send_mai_text('done with inpainting timing')
if 0 
	figure; plot(cumsum(time), err)
	hold on; plot(cumsum(time_circ), err_circ, 'r')
	figure; subplot(1,2,1); im(x); subplot(1,2,2); im(x - xtrue);
	figure; subplot(1,2,1); im(x_circ); subplot(1,2,2); im(x_circ - xtrue);
end

