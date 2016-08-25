% test tridiag inpaint
if ~isvar('wavelets')
	wavelets = 1;
end
fnames = {'autumn.tif'; 'bag.png'; 'cell.tif'; 'circuit.tif'; 'forest.tif'; ...
        'glass.png'; 'hands1.jpg'; 'lighthouse.png'; 'office_5.jpg'; ...
        'onion.png'; 'pout.tif'; 'rice.png'; 'toysnoflash.png'; 'trees.tif'};
comp = 0;
betas = [0.02 0.02 0.02 0.02 0.2 0.02 0.02 0.02 0.02 0.2 0.02 0.02 0.02 0.02];
sweep_range = 5;
for ii = 1%[4 6 10]
        fname = fnames{ii};
        xtrue = imread(fname);

        if size(xtrue, 3) > 1
                xtrue = rgb2gray(xtrue);
        end
        xtrue = double(xtrue);
        xtrue = downsample2(xtrue, 7);
        
        % make even dimensions
        if mod(size(xtrue, 1), 2) == 1
                xtrue = xtrue(2:end, :);
        end
        if mod(size(xtrue, 2), 2) == 1
                xtrue = xtrue(:, 2:end);
        end
        xtrue = xtrue(:,1:540);

        xtrue = xtrue./max(xtrue(:));
        
        [Nx, Ny] = size(xtrue);

	rng(0);
        reduce = 4;%8;
        samp = (rand(Nx, Ny) <= 1/reduce);
        D = Ginpaint(samp);
        [CH, CV] = construct_finite_diff([Nx Ny]);
	SNR = 40;
        y_noiseless = D * xtrue;
	sig = 10^(-SNR/20) * norm(y_noiseless) / sqrt(length(y_noiseless));
	y = y_noiseless + sig*randn(size(y_noiseless));
        %%
        betas = logspace(log10(reduce/100), log10(reduce/50), sweep_range);
        beta_big = 5;
%         beta = beta*(~samp) + beta_big*samp;
	% instead of zero-fill, let's do nearest-neighbor interp
	[xx, yy] = ndgrid(1:Nx, 1:Ny);
	xx_D = xx(samp);
	yy_D = yy(samp);
	xinit = griddata(xx_D, yy_D, y, xx, yy, 'nearest');
        niters = 200;
        mu0 = 1;%samp + 0.1;
        mu1 = mu0;
        
        R = [CH; CV];
        Rcirc = Cdiffs([Nx Ny],'offsets', [2 Nx], 'type_diff','circshift');
	if wavelets
        	betaws = logspace(log10(reduce/2000), log10(reduce/1000), sweep_range);
		W = Godwt1(true(Nx, Ny));
		%W = Gdiag(ones(Nx, Ny));
		RcircW = [Rcirc; W];
		
	else 
		betaws = 0;
		RcircW = Rcirc;
	end
        
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
			[x(:,:,jj,kk), xsaved, err(:,ii,jj,kk), cost, time(:,ii,jj,kk)] = AL_tridiag_inpaint(y, D, CHW, CVW, ...
					beta, xinit, xtrue, niters, 'mu', {1, 1, 1}, 'betaw', betaw, 'alphw', alphw);
if comp
			[x_P2(:,:,jj,kk), xsave_P2, err_P2(:,ii,jj,kk), costOrig_P2, time_P2(:,ii,jj,kk)] = AL_P2_inpainting(y, D, RW, ...
				xinit, niters, beta, xtrue, 'mu',  {1, 1} );
			
			[x_circ(:,:,jj), xsave_circ, err_circ(:,ii,jj,kk), costOrig_circ, time_circ(:,ii,jj,kk)] = AL_P2_inpainting(y, D, RcircW, ...
				xinit, niters, beta, xtrue, 'mu', {1, 1});
%   
end
		end
        end
%        best_beta_ndx = find(err_P2(end,ii,:) == min(err_P2(end,ii,:)))
        [best_beta_ndx, best_betaw_ndx] = find(squeeze(err(end,ii,:,:)) == min(col(err(end,ii,:,:))))
       if 0 
        [best_circ_beta_ndx, best_circ_betaw_ndx] = find(squeeze(err_circ(end,ii,:,:)) == min(col(err_circ(end,ii,:,:))))
		figure; plot(cumsum(time(:,ii,best_beta_ndx, best_betaw_ndx)), err(:,ii,best_beta_ndx, best_betaw_ndx))
		hold on; plot(cumsum(time_circ(:,ii,best_circ_beta_ndx, best_circ_betaw_ndx)), err_circ(:,ii,best_circ_beta_ndx, best_circ_betaw_ndx), 'r')
		%         keyboard
		figure; subplot(1,2,1); im(x(:,:,best_beta_ndx)); subplot(1,2,2); im(x(:,:,best_beta_ndx) - xtrue);
		figure; subplot(1,2,1); im(x_circ(:,:,best_circ_beta_ndx)); subplot(1,2,2); im(x_circ(:,:,best_circ_beta_ndx) - xtrue);
	end
end

save(sprintf('inpainting_%s_wavelet%d', machine(1:3), wavelets))
