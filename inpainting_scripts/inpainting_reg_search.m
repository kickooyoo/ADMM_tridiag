% test tridiag inpaint
sweep_range = 7 - 2*wavelets;
betas = logspace(log10(reduce/500), log10(reduce/5), sweep_range);
if wavelets
betaws = logspace(log10(reduce/200), log10(reduce/50), sweep_range);
else
	betaws = 0;
end
betas_circ = betas;
betaws_circ = betaws;
min_beta_factor = 1e-6;
search_range = 100;
save_fname = sprintf('inpainting_mat/inpainting_reg_%s_wavelet%d_SNR%d_reduce%1.2d.mat', machine(1:3), wavelets, SNR, reduce);
err_mask = true(Nx, Ny);
mask_width = 1;
err_mask([1:mask_width, end-mask_width+1:end],:) = false;
err_mask(:,[1:mask_width, end-mask_width+1:end]) = false;

for ii = 1:10
	for jj = 1:length(betas)
		beta = betas(jj);
		for kk = 1:length(betaws)
			if wavelets
				betaw = betaws(kk);
				alphw = 1;
				CHW = [CH; betaw * alphw / beta * W];
				CVW = [CV; betaw * (1-alphw) / beta * W];     
				RW = [CHW; CVW];
			else
				betaw = 0;
				CHW = CH;
				CVW = CV;
				RW = R;
			end
			[x(:,:,jj,kk), xsaved, err(:,jj,kk), cost, time(:,jj,kk)] = AL_tridiag_inpaint(y, D, CHW, CVW, ...
					beta, xinit, xtrue, niters, 'mu', {mu0, mu1, mu2}, 'betaw', betaw, 'alphw', alphw, 'mask', err_mask);
		end
	end
	save(save_fname);

	% adjust betas, betaws
	if wavelets
		[best_beta_ndx, best_betaw_ndx] = find(squeeze(err(end,:,:)) == min(col(err(end,:,:))), 1);
	else
		best_beta_ndx = find(squeeze(err(end,:,:)) == min(col(err(end,:,:))), 1);
	end
	if best_beta_ndx == 1
		left_beta_factor = max(reduce./betas(1) - search_range, min_beta_factor);
	else
		left_beta_factor = reduce./betas(best_beta_ndx - 1);
	end
	if best_beta_ndx == length(betas)
		right_beta_factor = reduce./betas(end) + search_range;
	else
		right_beta_factor = reduce./betas(best_beta_ndx + 1);
	end
	betas = logspace(log10(reduce/left_beta_factor), log10(reduce/right_beta_factor), sweep_range);
	if wavelets
		if best_betaw_ndx == 1
			left_betaw_factor = reduce./max(betaws(1) - search_range, min_beta_factor);
		else
			left_betaw_factor = reduce./betaws(best_betaw_ndx - 1);
		end
		if best_betaw_ndx == length(betaws)
			right_betaw_factor = reduce./(betaws(end) + search_range);
		else
			right_betaw_factor = reduce./betaws(best_betaw_ndx + 1);
		end
		betaws = logspace(log10(reduce/left_betaw_factor), log10(reduce/right_betaw_factor), sweep_range);
	end
	for jj = 1:length(betas_circ)
		beta_circ = betas_circ(jj);
		for kk = 1:length(betaws_circ)
			if wavelets
				betaw = betaws_circ(kk);
				alphw = 1;
				CHW = [CH; betaw * alphw / beta_circ * W];
				CVW = [CV; betaw * (1-alphw) / beta_circ * W];     
				RW = [CHW; CVW];
			else
				betaw = 0;
				CHW = CH;
				CVW = CV;
				RW = R;
			end
			[x_circ(:,:,jj,kk), xsave_circ, err_circ(:,jj,kk), costOrig_circ, time_circ(:,jj,kk)] = AL_P2_inpainting(y, D, RcircW, ...
				xinit, niters, beta_circ, xtrue, 'mu', {mu0, mu1}, 'mask', err_mask);
		end
	end
	save(save_fname);

	% adjust betas_circ, betaws_circ
	if wavelets
		[best_circ_beta_ndx, best_circ_betaw_ndx] = find(squeeze(err_circ(end,:,:)) == min(col(err_circ(end,:,:))), 1);
	else
		best_circ_beta_ndx = find(squeeze(err_circ(end,:,:)) == min(col(err_circ(end,:,:))), 1);
	end
	if best_circ_beta_ndx == 1
		left_beta_factor = max(reduce./betas_circ(1) - search_range, min_beta_factor);
	else
		left_beta_factor = reduce./betas_circ(best_circ_beta_ndx - 1);
	end
	if best_circ_beta_ndx == length(betas_circ)
		right_beta_factor = reduce./betas_circ(end) + search_range;
	else
		right_beta_factor = reduce./betas_circ(best_circ_beta_ndx + 1);
	end
	betas_circ = logspace(log10(reduce/left_beta_factor), log10(reduce/right_beta_factor), sweep_range);
	if wavelets
		if best_circ_betaw_ndx == 1
			left_betaw_factor = reduce./max(betaws(1) - search_range, min_beta_factor);
		else
			left_betaw_factor = reduce./betaws(best_circ_betaw_ndx - 1);
		end
		if best_circ_betaw_ndx == length(betaws)
			right_betaw_factor = reduce./(betaws(end) + search_range);
		else
			right_betaw_factor = reduce./betaws(best_circ_betaw_ndx + 1);
		end
		betaws_circ = logspace(log10(reduce/left_betaw_factor), log10(reduce/right_betaw_factor), sweep_range);
	end
	display(sprintf('done with %d rounds of coarse to fine', ii));
end


