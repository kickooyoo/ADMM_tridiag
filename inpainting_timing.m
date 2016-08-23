% test tridiag inpaint
inpainting_setup;

sweep_range = 5;
betas = logspace(log10(reduce/100), log10(reduce/50), sweep_range);
betaws = logspace(log10(reduce/2000), log10(reduce/1000), sweep_range);

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
				beta, xinit, xtrue, niters, 'mu', {mu0, mu1, mu2}, 'betaw', betaw, 'alphw', alphw);
		[x_P2(:,:,jj,kk), xsave_P2, err_P2(:,jj,kk), costOrig_P2, time_P2(:,jj,kk)] = AL_P2_inpainting(y, D, RW, ...
			xinit, niters, beta, xtrue, 'mu',  {mu0, mu1} );
		
		[x_circ(:,:,jj), xsave_circ, err_circ(:,jj,kk), costOrig_circ, time_circ(:,jj,kk)] = AL_P2_inpainting(y, D, RcircW, ...
			xinit, niters, beta, xtrue, 'mu', {mu0, mu1});
%   
	end
end
save(sprintf('inpainting_timing_%s_wavelet%d', machine(1:3), wavelets))

[best_beta_ndx, best_betaw_ndx] = find(squeeze(err(end,:,:)) == min(col(err(end,:,:))))
[best_circ_beta_ndx, best_circ_betaw_ndx] = find(squeeze(err_circ(end,:,:)) == min(col(err_circ(end,:,:))))
if 0 
	figure; plot(cumsum(time(:,best_beta_ndx, best_betaw_ndx)), err(:,best_beta_ndx, best_betaw_ndx))
	hold on; plot(cumsum(time_circ(:,best_circ_beta_ndx, best_circ_betaw_ndx)), err_circ(:,best_circ_beta_ndx, best_circ_betaw_ndx), 'r')
	%         keyboard
	figure; subplot(1,2,1); im(x(:,:,best_beta_ndx)); subplot(1,2,2); im(x(:,:,best_beta_ndx) - xtrue);
	figure; subplot(1,2,1); im(x_circ(:,:,best_circ_beta_ndx)); subplot(1,2,2); im(x_circ(:,:,best_circ_beta_ndx) - xtrue);
end

