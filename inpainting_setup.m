% test tridiag inpaint

fnames = {'autumn.tif'; 'bag.png'; 'cell.tif'; 'circuit.tif'; 'forest.tif'; ...
        'glass.png'; 'hands1.jpg'; 'lighthouse.png'; 'office_5.jpg'; ...
        'onion.png'; 'pout.tif'; 'rice.png'; 'toysnoflash.png'; 'trees.tif'};
betas = [0.02 0.02 0.02 0.02 0.2 0.02 0.02 0.02 0.02 0.2 0.02 0.02 0.02 0.02];
for ii = 6%[4 6 10]
        fname = fnames{ii};
        xtrue = imread(fname);
        if size(xtrue, 3) > 1
                xtrue = rgb2gray(xtrue);
        end
        % make even dimensions
        if mod(size(xtrue, 1), 2) == 1
                xtrue = xtrue(2:end, :);
        end
        if mod(size(xtrue, 2), 2) == 1
                xtrue = xtrue(:, 2:end);
        end
        xtrue = double(xtrue);
        xtrue = xtrue./max(xtrue(:));
        [Nx, Ny] = size(xtrue);
        rng(0);
        reduce = 4;%8;
        samp = (rand(Nx, Ny) <= 1/reduce);
        D = Ginpaint(samp);
        [CH, CV] = construct_finite_diff([Nx Ny]);
        y = D * xtrue;
        %%
        beta_smalls = logspace(log10(reduce/350), log10(reduce/250), 10);
        beta_big = 5;
%         beta = beta_small*(~samp) + beta_big*samp;
        xinit = reshape(D'*y, Nx, Ny);
        niters = 200;
        mu0 = 1;%samp + 0.1;
        mu1 = mu0;
        %%
        % [x, xsaved, err, cost, time] = AL_tridiag_inpaint(y, D, CH, CV, ...
        %         beta, xinit, xtrue, niters, 'mu', {mu0, mu1, 1});
        % figure; im(x)
        % calc_NRMSE_over_mask(x, xtrue, true(size(x)))
        % calc_NRMSE_over_mask(x, xtrue, ~samp)
        %
        % [x2, xsaved, err, cost, time] = AL_tridiag_inpaint(y, D, CH, CV, ...
        %         beta_small, xinit, xtrue, niters, 'mu', {mu0, mu1, 1});
        % figure; im(x2)
        % calc_NRMSE_over_mask(x2, xtrue, true(size(x)))
        % calc_NRMSE_over_mask(x2, xtrue, ~samp)
        %
        % beta = 10;
        % [x3, xsaved, err, cost, time] = AL_tridiag_inpaint(y, D, CH, CV, ...
        %         beta_big, xinit, xtrue, niters, 'mu', {mu0, mu1, 1});
        % figure; im(x3)
        % calc_NRMSE_over_mask(x3, xtrue, true(size(x)))
        % calc_NRMSE_over_mask(x2, xtrue, ~samp)
        
        %%
        R = [CH; CV];
        Rcirc = Cdiffs([Nx Ny],'offsets', [1 Nx], 'type_diff','circshift');
        
        for jj = 1:length(beta_smalls)
                beta_small = beta_smalls(jj);
        [x(:,:,jj), xsave, err(:,ii,jj), costOrig, time(:,ii,jj)] = AL_P2_inpainting(y, D, R, ...
                xinit, niters, beta_small, xtrue, 'mu', {1, 1});
%         figure; subplot(1,2,1); im(x); subplot(1,2,2); im(x - xtrue);
        
        [x_P2, xsave_P2, err_P2(:,ii,jj), costOrig_P2, time_P2(:,ii,jj)] = AL_P2_inpainting(y, D, R, ...
                xinit, niters, beta_small, xtrue, 'mu',  {1, 1});
        

        [x_circ(:,:,jj), xsave_circ, err_circ(:,ii,jj), costOrig_circ, time_circ(:,ii,jj)] = AL_P2_inpainting(y, D, Rcirc, ...
                xinit, niters, beta_small, xtrue, 'mu', {1, 1});
%         

        end
        best_beta_ndx = find(err_P2(end,ii,:) == min(err_P2(end,ii,:)))
        best_beta_ndx = find(err(end,ii,:) == min(err(end,ii,:)))
        best_circ_beta_ndx = find(err_circ(end,ii,:) == min(err_circ(end,ii,:)))
        
        figure; plot(cumsum(time(:,ii,best_beta_ndx)), err(:,ii,best_beta_ndx))
        hold on; plot(cumsum(time_circ(:,ii,best_circ_beta_ndx)), err_circ(:,ii,best_circ_beta_ndx), 'r')
        %         keyboard
        figure; subplot(1,2,1); im(x(:,:,best_beta_ndx)); subplot(1,2,2); im(x(:,:,best_beta_ndx) - xtrue);
        figure; subplot(1,2,1); im(x_circ(:,:,best_circ_beta_ndx)); subplot(1,2,2); im(x_circ(:,:,best_circ_beta_ndx) - xtrue);
end


