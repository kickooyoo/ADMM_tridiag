exp_setup;


%method = 'tridiag';
method = 'MFISTA';
if gen
	niters = 7500;
	% already done 2.^(3:20);
	betas = 2.^(19:0.5:23);
	betas = 2.^(18:27);
%	[bc_sc, sc] = ir_wls_init_scale(F*S, y_noise, body_coil);
	xinit_tmp = xinit;%zeros(size(xinit));%(xinit + bc_sc)./2;
else
	betas = 2.^(16:31);
end
xhat_betas = [];
for ii = 1:length(betas)
	beta = betas(ii);

	switch method
		case 'tridiag'
			fname = sprintf('%s/x_tri_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta);
			if gen
				[x_tri_inf, ~, ~, costOrig_tri, time_tri] = ADMM_tridiag(y_noise, F, S, CH, CV, beta, xinit_tmp, zeros(size(SoS)), niters);
				save(fname, 'x_tri_inf', 'body_coil');
				xhat_betas(:,:,ii) = x_tri_inf;
				body_coil_err(ii) = calc_NRMSE(x_tri_inf./max(abs(col(x_tri_inf))), body_coil./max(abs(col(body_coil))), mask);
			elseif exist(fname)
				load(fname, 'x_tri_inf', 'body_coil');
				xhat_betas(:,:,ii) = x_tri_inf;
				body_coil_err(ii) = calc_NRMSE(x_tri_inf./max(abs(col(x_tri_inf))), body_coil./max(abs(col(body_coil))), mask);
			else
				display(sprintf('%s does not exist (beta = 2^%d)', fname, log2(beta)))
			end
		case 'MFISTA' 
			fname = sprintf('%s/x_MFISTA_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta);
			if gen
				x_MFISTA = MFISTA_wrapper(Nx, Ny, R, y_noise, xinit_tmp, F, S, beta, niters, curr_folder);
				save(fname, 'x_MFISTA', 'body_coil');
				%if isvar('x_MFISTA')
					xhat_betas(:,:,ii) = x_MFISTA;
				%elseif isvar('xMFIS')
				%	xhat_betas(:,:,ii) = xMFIS;
				%else 
				%	display('where is MFISTA output?')
				%end
				body_coil_err(ii) = calc_NRMSE(x_MFISTA./max(abs(col(x_MFISTA))), body_coil./max(abs(col(body_coil))), mask);
			elseif exist(fname)
				load(fname, '*MFIS*', 'body_coil');
				%if isvar('x_MFISTA')
					xhat_betas(:,:,ii) = x_MFISTA;
				%elseif isvar('xMFIS')
				%	xhat_betas(:,:,ii) = xMFIS;
				%else 
				%	display('where is MFISTA output?')
				%end
				body_coil_err(ii) = calc_NRMSE(x_MFISTA./max(abs(col(x_MFISTA))), body_coil./max(abs(col(body_coil))), mask);
			else
				display(sprintf('%s does not exist (beta = 2^%d)', fname, log2(beta)))
			end
		otherwise
			keyboard;
	end
end
figure; im(xhat_betas);
if gen, send_mai_text(sprintf('done searching betas on %s, time to choose beta and gen xinf', machine(1:3))); end
