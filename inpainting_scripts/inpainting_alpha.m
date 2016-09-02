inpainting_setup;

alphas = [0 0.25 0.5 0.75 1];
wavelets = 0;
niters = 2000;

for jj = 1:length(alphas)
	alph = alphas(jj);
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
				beta, xinit, xtrue, niters, 'mu', {1, 1, 1}, 'betaw', betaw, 'alphw', alphw, 'alph', alph);

	end
end

save(sprintf('inpainting_alpha_%s_wavelet%d', machine(1:3), wavelets))
