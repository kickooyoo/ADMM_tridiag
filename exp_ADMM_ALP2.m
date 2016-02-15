% experimental data for tridiag_ADMM

tridiag_exp_setup;
niters = 350;

if 1
	% load xinfs
	load(sprintf('./reviv/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_tri_inf');
	if truncate
		x_tri_inf = reshape(x_tri_inf, 256, 144);
		x_tri_inf = x_tri_inf(3:end-2, 3:end-2);
	else
		x_tri_inf = x_tri_inf;
	end
	display('using tridiag ADMM as true!')
else 
	load(sprintf('./reviv/x_MFISTA_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x*MFIS*');
		if isvar('xMFIS')
			x_MFISTA = xMFIS;
		end
	if truncate
		xMFIS = reshape(x_MFISTA, 256, 144);
		x_tri_inf = xMFIS(3:end-2, 3:end-2);
	else
		x_tri_inf = x_MFISTA;
	end
end

% maybe not xinf since not guaranteed to converge
% load([home 'Documents/mai_code/ADMM_tridiag/reviv/al_p2_chcv_5000iter.mat'],'xhat_alp2');
% if truncate
%         xhat_alp2 = reshape(xhat_alp2, 256, 144);
%         x_alp2_inf = xhat_alp2(3:end-2, 3:end-2);
% else
%         x_alp2_inf = xhat_alp2;
% end

if 0
load(sprintf('./reviv/x_alp2c_inf_slice%d_beta%.*d.mat', slice, 3, beta),'x_alp2c_inf');
if truncate
        x_alp2c_inf = reshape(x_alp2c_inf, 256, 144);
        x_alp2c_inf = x_alp2c_inf(3:end-2, 3:end-2);
else
        x_alp2c_inf = x_alp2c_inf;
end
end

if wavelets 
[xhat_alp2c, ~, nrmsd_alp2c, costOrig_alp2c, time_alp2c] = AL_P2_gen(...
        y_noise, F, S, RcircW, SoS, niters, beta, x_alp2c_inf, 'zmethod', 'fft');

[xhat_alp2t, ~, nrmsd_alp2t, costOrig_alp2t,time_alp2t] = AL_P2_gen(...
        y_noise, F, S, RW, SoS, niters, beta, x_tri_inf, 'inner_iter', 3);
end

for ii = 1:1%3%5
        [xhat_alp2t(:,:,ii), ~, nrmsd_alp2t(:,ii), costOrig_alp2t(:,ii), time_alp2t(:,ii)] = AL_P2_gen(y_noise, F, S, R, xinit, niters, beta, x_tri_inf,'inner_iter', ii, 'mask', mask);
	if (length(time_alp2t) ~= niters + 1), keyboard; end
end

%[xhat_alp2c, ~, nrmsd_alp2c, costOrig_alp2c, time_alp2c] = AL_P2_gen(y_noise, F, S, Rcirc, xinit, niters, beta, x_alp2c_inf, 'zmethod','fft', 'mask', mask);
[xhat_tri, ~, nrmsd_tri, costOrig_tri, time_tri] = tridiag_ADMM(y_noise, F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'mu', ones(1,5));
%[xhat_tri_max, ~, nrmsd_tri_max, costOrig_tri_max, time_tri_max] = tridiag_ADMM(y_noise, F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'nthread', int32(maxNumCompThreads('automatic')),'mu', ones(1,5));
%[xhat_tri_2max, ~, nrmsd_tri_2max, costOrig_tri_2max, time_tri_2max] = tridiag_ADMM(y_noise, F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'nthread', int32(160),'mu', ones(1,5));
%[xhat_tri_mu, ~, nrmsd_tri_mu, costOrig_tri_mu, time_tri_mu] = tridiag_ADMM(y_noise, F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask);
%[xhat_tri_max_mu, ~, nrmsd_tri_max_mu, costOrig_tri_max_mu, time_tri_max_mu] = tridiag_ADMM(y_noise, F, S, CH, CV, beta, xinit, x_tri_inf, niters, 'mask', mask, 'nthread', int32(maxNumCompThreads('automatic')));

save(sprintf('./reviv/mpel8_timing_%dx%d_%diter_%dslice_ADMMtrue.mat', Nx, Ny, niters, slice));
send_mai_text('done with mpel8 timing');

display('DONE');
