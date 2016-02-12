tridiag_exp_setup;
niters = 10;
x = xinit;

rng(0);
num_test = 100;

rand_scale = round(100*rand(1, 1, num_test));
rand_noise = 100*randn(Nx, Ny, num_test);
test_xs = repmat(rand_scale, [Nx Ny 1]).*repmat(x,[1 1 num_test]) + rand_noise; 

test_y = 0;
betas = 2.^(1:30);
test_F = 0;
test_S = 0;
for jj = 1:length(betas)
	beta_test = betas(jj);
	%for ii = 1:num_test
		x_test = test_xs(:,:,ii);
		cost_tridiag(ii,jj) = tridiag_ADMM_cost(test_y, test_F, test_S, CH, CV, beta_test, x_test);

		cost_MFISTA(ii,jj) = MFISTA_cost_wrapper(Nx, Ny, R, test_y, x_test, test_F, test_S, beta_test);
	%end
end

figure; plot(col(cost_tridiag)./col(cost_MFISTA))
figure; plot(betas, cost_tridiag, 'r'); hold on; plot(betas, cost_MFISTA,'b');
figure; scatter(cost_tridiag, cost_MFISTA)



