if gen
	tridiag_exp_setup;
else
	slice = 67;
	slice = 38;
end
niters = 1000;

% already done 2.^(3:20);
betas = 2.^(20:30);

method = 'tridiag';
method = 'MFISTA';

clear xhat_betas
for ii = 1:length(betas)
	beta = betas(ii);

	if gen 
		switch method
			case 'tridiag'
				[xhat_betas(:,:,ii), ~, ~, costOrig_tri, time_tri] = tridiag_ADMM(y, F, S, CH, CV, alph, beta, xinit, zeros(size(SoS)), niters, 'mu', mu);
				x_tri_inf = xhat_tri;
				save(sprintf('reviv/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_tri_inf');
			case 'MFISTA' 
				x_MFISTA = MFISTA_wrapper(Nx, Ny, R, y, xinit, F, S, beta, niters);
				xhat_betas(:,:,ii) = x_MFISTA;
				save(sprintf('reviv/x_MFISTA_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_MFISTA');
			otherwise
				keyboard;
		end
	else
		switch method
			case 'tridiag'
				fname = sprintf('reviv/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta);
				if exist(fname)
					load(fname, 'x_tri_inf');
					xhat_betas(:,:,ii) = x_tri_inf;
				else
					display(sprintf('%s does not exist (beta = 2^%d)', fname, log2(beta)))
				end
			case 'MFISTA'
				fname = sprintf('reviv/x_MFISTA_inf_slice%d_beta%.*d.mat', slice, 3, beta);
				if exist(fname)
					load(fname, 'x_MFISTA');
					xhat_betas(:,:,ii) = x_MFISTA;
				else
					display(sprintf('%s does not exist (beta = 2^%d)', fname, log2(beta)))
				end
			otherwise
				keyboard;
		end
	end
end
figure; jf_slicer(xhat_betas);
figure; im(xhat_betas);
if gen, send_mai_text('done searching betas'); end
