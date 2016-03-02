if gen
	exp_setup;
else
	switch orient
		case 'axial'
			slice = 38;
		case 'sagittal'
			slice = 69;
		case 'coronal'
			slice = 135; 
		otherwise
			display(sprintf('need to choose slice for %s orientation', orient))
			keyboard
	end
end
niters = 700;

% already done 2.^(3:20);
betas = 2.^(19:0.5:23);

%method = 'tridiag';
method = 'MFISTA';
%if ~issim && ~isvar('body_coil')
%	load(sprintf('body_coil_slice%d.mat', slice),'body_coil')
%end

%if issim
%	[sense_maps, body_coil, Sxtrue] = sim_setup();
%	[Nx, Ny, Nc] = size(Sxtrue);
%	mask = true(Nx, Ny);
%end 

clear xhat_betas
for ii = 1:length(betas)
	beta = betas(ii);

	switch method
		case 'tridiag'
			fname = sprintf('%s/x_tri_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta);
			if gen
				[x_tri_inf, ~, ~, costOrig_tri, time_tri] = ADMM_tridiag(y, F, S, CH, CV, beta, xinit, zeros(size(SoS)), niters);
			elseif exist(fname)
				load(fname, 'x_tri_inf', 'body_coil');
			else
				display(sprintf('%s does not exist (beta = 2^%d)', fname, log2(beta)))
			end
			xhat_betas(:,:,ii) = x_tri_inf;
			curr_img = xhat_betas(:,:,ii);
			body_coil_err(ii) = calc_NRMSE_over_mask(curr_img./max(abs(col(curr_img))), body_coil./max(abs(col(body_coil))), mask);
			if gen
				save(fname, 'x_tri_inf', 'body_coil');
			end
		case 'MFISTA' 
			fname = sprintf('%s/x_MFISTA_inf_%s_beta%.*d.mat', curr_folder, slice_str, 3, beta);
			if gen
				x_MFISTA = MFISTA_wrapper(Nx, Ny, R, y, xinit, F, S, beta, niters);
			elseif exist(fname)
				load(fname, 'x_MFISTA', 'body_coil');
			else
				display(sprintf('%s does not exist (beta = 2^%d)', fname, log2(beta)))
			end
			xhat_betas(:,:,ii) = x_MFISTA;
			curr_img = xhat_betas(:,:,ii);
			body_coil_err(ii) = calc_NRMSE_over_mask(curr_img./max(abs(col(curr_img))), body_coil./max(abs(col(body_coil))), mask);
			if gen
				save(fname, 'x_MFISTA', 'body_coil');
			end
		otherwise
			keyboard;
	end
end
figure; im(xhat_betas);
if gen, send_mai_text('done searching betas, time to choose beta and gen xinf'); end
