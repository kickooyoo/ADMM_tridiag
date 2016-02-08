if gen
	tridiag_exp_setup;
else
	slice = 67;
	slice = 38;
end
niters = 1000;

% already done 2.^(3:20);
betas = 2.^(3:30);

for ii = 1:length(betas)
	beta = betas(ii);

	if gen 
		%[xhat_betas(:,:,ii), ~, ~, costOrig_tri, time_tri] = tridiag_ADMM(y, F, S, CH, CV, alph, beta, SoS, zeros(size(SoS)), niters, 'mu', mu);
		%x_tri_inf = xhat_tri;
		xhat_betas(:,:,ii) = MFISTA_wrapper(Nx, Ny, R, y, xinit, F, S, beta, niters)
		save(sprintf('./reviv/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_tri_inf');
	else
		load(sprintf('./reviv/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_tri_inf');
		xhat_betas(:,:,ii) = x_tri_inf;
	end
end
if gen
	figure; jf_slicer(xhat_betas);
end
send_mai_text('done searching betas')
