% test tridiag inpaint
if ~isvar('niters')
	niters = 20000;
end
if ~isvar('do_alph')
	do_alph = 0;
end
if ~isvar('str_mod')
	str_mod = '';
end
if do_alph
	alphas = 0:0.1:1;
	save_fname = sprintf('inpainting_mat/%s/timing/inpainting_timing_%s_iters%d_wavelet%d_SNR%d_reduce%1.2d_%strue_alphw%1.1d_retunedmu%s.mat', obj, machine(1:3), niters, wavelets, SNR, reduce, true_opt, alphw, str_mod);
else
	alphas = 0.5;
	save_fname = sprintf('inpainting_mat/%s/timing/inpainting_timing_%s_iters%d_wavelet%d_SNR%d_reduce%1.2d_%strue_%1.1dalph_tunedmu%s.mat', obj, machine(1:3), niters, wavelets, SNR, reduce, true_opt, alphas, str_mod);
end
if exist(save_fname, 'file')
	display('file already exists');
	keyboard;
end

do_exp = [0 1 0 1 1 0];
if do_exp(1)
	aaa = find(alphas == 0.5);
	[x(:,:,aaa), ~, err(:,aaa), cost(:,aaa), time(:,aaa)] = AL_tridiag_inpaint(y, D, CHW, CVW, ...
			beta, xinit, xtrue, niters, 'betaw', betaw, 'alphw', alphw, 'alph', alph, 'kapp', 7, 'pos', 0.01);
end
if do_exp(2)
	[x_P2, ~, err_P2, costOrig_P2, time_P2] = AL_P2_inpainting(y, D, RW, ...
		xinit, niters, beta, xtrue, 'alphw', alphw, 'betaw', betaw, 'inner_iter', 5, 'save_progress', [save_fname '_tmp']);
end
if do_exp(3)
	[x_circ, ~, err_circ, costOrig_circ, time_circ] = AL_P2_inpainting(y, D, RcircW, ...
		xinit, niters, beta, xtrue, 'alphw', alphw, 'betaw', betaw, 'zmethod', 'FFT');
end
if do_exp(4)
	[x_MFIS, C_MFIS, time_MFIS, err_MFIS, ~] = MFISTA_inpainting_wrapper(Nx, Ny, RW, y, xinit, D, beta, round(niters/2), curr_folder, slice_str, 'xinf', xtrue, 'xinfnorm', norm(col(xtrue),2));
	save(save_fname, '-append')
end
if do_exp(5)
	[x_ADMM, ~, err_ADMM, cost_ADMM, time_ADMM] = ADMM_tridiag_inpaint(y, D, CHW, CVW, beta, xinit, xtrue, niters, 'betaw', betaw, 'alphw', alphw, 'alph', alph, 'kapp', 7, 'pos', 0.01);
end
if exist(save_fname, 'file')
	save(save_fname, '-append')
else
	save(save_fname)
end
send_mai_text('done with inpainting timing except alph')
if do_exp(6)
	for aa = setdiff(1:length(alphas), aaa)
		alph = alphas(aa);
		[x(:,:,aa), ~, err(:,aa), cost(:,aa), time(:,aa)] = AL_tridiag_inpaint(y, D, CHW, CVW, ...
			beta, xinit, xtrue, niters, 'betaw', betaw, 'alphw', alphw, 'alph', alph, 'kapp', 7, 'pos', 0.01);
	end
save(save_fname,'-append')
end
send_mai_text('done with inpainting timing')

