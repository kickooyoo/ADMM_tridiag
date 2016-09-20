
alphaws = [0 0.25 0.5 0.75 1];
wavelets = 1;
niters = 2000;

if ~isvar('str_mod')
	str_mod = '';
end

save_fname = sprintf('inpainting_mat/%s/timing/inpainting_alphaw_%s_iters%d_wavelet%d_SNR%d_reduce%1.2d_%strue_retunedmu%s.mat', obj, machine(1:3), niters, wavelets, SNR, reduce, true_opt, str_mod);
for jj = 1:length(alphaws)
	alphw = alphaws(jj);
	if alphw == 0
		CHW = CH; 
		CVW = [CV; betaw / beta * W]; 
	elseif alphw == 1
		CHW = [CH; betaw / beta * W]; 
		CVW = CV; 
	else
		CHW = [CH; betaw * alphw / beta * W]; 
		CVW = [CV; betaw * (1-alphw) / beta * W]; 
	end 
	%RW = [CH; CV; betaw / beta * W]; 
	%RcircW = [Rcirc; betaw_circ / beta_circ * W]; 
	[x(:,:,jj), ~, err(:,jj), cost, time(:,jj)] = AL_tridiag_inpaint(y, D, CHW, CVW, ...
		beta, xinit, xtrue, niters, 'mu', {1, 1, 1}, 'betaw', betaw, 'alphw', alphw, 'alph', alph);
	display(sprintf('done with alphaw = %d', alphw))	
end

