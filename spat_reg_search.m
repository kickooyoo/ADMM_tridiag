if gen
	exp_setup;
else
	if ~issim
		slice = 67;
		slice = 38;
	end
end
niters = 1000;

% already done 2.^(3:20);
betas = 2.^(13:21);

method = 'tridiag';
%method = 'MFISTA';
if ~issim && ~isvar('body_coil')
	load(sprintf('body_coil_slice%d.mat', slice),'body_coil')
end

if issim
	[sense_maps, body_coil, Sxtrue] = sim_setup();
	[Nx, Ny, Nc] = size(Sxtrue);
	mask = true(Nx, Ny);
end
if issim 
	slice_str = 'sim';
else
	slice_str = sprintf('slice%d', slice);
end

clear xhat_betas
for ii = 1:length(betas)
	beta = betas(ii);

	if gen 
		switch method
			case 'tridiag'
				[x_tri_inf, ~, ~, costOrig_tri, time_tri] = ADMM_tridiag(y, F, S, CH, CV, beta, xinit, zeros(size(SoS)), niters);
				xhat_betas(:,:,ii) = x_tri_inf;
				save(sprintf('./reviv/curr/x_tri_inf_%s_beta%.*d.mat', slice_str, 3, beta), 'x_tri_inf');
			case 'MFISTA' 
				x_MFISTA = MFISTA_wrapper(Nx, Ny, R, y, xinit, F, S, beta, niters);
				xhat_betas(:,:,ii) = x_MFISTA;
				save(sprintf('./reviv/curr/x_MFISTA_inf_%s_beta%.*d.mat', slice_str, 3, beta), 'x_MFISTA');
			otherwise
				keyboard;
		end
	else
		switch method
			case 'tridiag'
				fname = sprintf('./reviv/curr/x_tri_inf_%s_beta%.*d.mat', slice_str, 3, beta);
				if exist(fname)
					load(fname, 'x_tri_inf');
					xhat_betas(:,:,ii) = x_tri_inf;
				else
					display(sprintf('%s does not exist (beta = 2^%d)', fname, log2(beta)))
				end
			case 'MFISTA'
				fname = sprintf('./reviv/curr/x_MFISTA_inf_%s_beta%.*d.mat', slice_str, 3, beta);
				if exist(fname)
					load(fname, 'x_MFISTA');
					xhat_betas(:,:,ii) = x_MFISTA;
				else
					display(sprintf('%s does not exist (beta = 2^%d)', fname, log2(beta)))
				end
			otherwise
				keyboard;
		end
		curr_img = xhat_betas(:,:,ii);
		body_coil_err(ii) = calc_NRMSE_over_mask(curr_img./max(abs(col(curr_img))), body_coil./max(abs(col(body_coil))), mask);
	end
end
figure; im(xhat_betas);
if gen, send_mai_text('done searching betas'); end
