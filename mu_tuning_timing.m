% experimental data for tridiag_ADMM

tridiag_exp_setup;
niters = 50;

% cost function for MFISTA lacks 1/2
load(sprintf('./reviv/x_MFISTA_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x*MFIS*');
if isvar('xMFIS')
        x_MFISTA = xMFIS;
end
x_tri_inf = x_MFISTA;

load(sprintf('./reviv/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_tri_inf');

[xhat_trit, ~, nrmsd_trit, costOrig_trit, time_trit] = tridiag_ADMM(y_noise, ...
        F, S, CH, CV, alph, beta, xinit, x_tri_inf, niters, 'mask', mask);
[xhat_trit2, ~, nrmsd_trit2, costOrig_trit2, time_trit2] = tridiag_ADMM(y_noise, ...
        F, S, CH, CV, alph, beta, xinit, x_tri_inf, niters, 'mask', mask, 'fancy_mu34', 1);

[xhat_triat2, ~, nrmsd_triat2, costOrig_triat2, time_triat2] = tridiag_ADMM(y_noise, ...
        F, S, CH, CV, 0, beta, xinit, x_tri_inf, niters, 'mask', mask, 'fancy_mu34', 1);

plain_mu = num2cell(ones(1,5));
[xhat_tri, ~, nrmsd_tri, costOrig_tri, time_tri] = tridiag_ADMM(y_noise, ...
        F, S, CH, CV, alph, beta, xinit, x_tri_inf, niters, 'mask', mask, 'mu', plain_mu);
figure; plot(nrmsd_trit); hold on; plot(nrmsd_trit2,'r'); 
hold on; plot(nrmsd_triat2,'k'); hold on; plot(nrmsd_tri,'m');

[xhat_alp2, ~, nrmsd_alp2, costOrig_alp2, time_alp2] = AL_P2_gen(y_noise, ...
        F, S, R, xinit, niters, beta, x_tri_inf,'inner_iter', 1, 'mask', mask, 'mu', plain_mu(1:3));

[xhat_alp2t, ~, nrmsd_alp2t, costOrig_alp2t, time_alp2t] = AL_P2_gen(y_noise, ...
        F, S, R, xinit, niters, beta, x_tri_inf,'inner_iter', 1, 'mask', mask);



save(sprintf('./reviv/mu_test_%dx%d_%diter_%dslice.mat', Nx, Ny, niters, slice));
send_mai_text('done with mu exp timing');

display('DONE');
figure; plot(cumsum(time_trit), nrmsd_trit); 
hold on; plot(cumsum(time_alp2t), nrmsd_alp2t,'r')
plot(cumsum(time_trit2), nrmsd_trit2, 'g');
plot(cumsum(time_tri), nrmsd_tri, 'm');

%%
return
%%
[xhat_alp2t, ~, nrmsd_alp2t, costOrig_alp2t, time_alp2t] = AL_P2_gen(y_noise, ...
        F, S, R, xinit, niters, beta, x_tri_inf,'inner_iter', 1, 'mask', mask, 'parFFT', true);

%%
slice = 38;
load(sprintf('%sDocuments/data/2010-07-06-fessler-3d/slice38/ramani/Smaps%d.mat', home_path, slice), 'Smap_QPWLS');
SS = sos_combine(permute(Smap_QPWLS, [1 2 4 3]));
maxSS = max(col(SS));
minSS = min(col(SS));

%%
Nx = 5; Ny = 7;
[CH, CV] = construct_finite_diff([Nx Ny]);
R = [CH; CV];
Rcirc = Cdiffs([Nx Ny],'offsets', [1 Nx], 'type_diff','circshift');
CH2 = full(CH'*CH);
CHcirc = Cdiffs([Nx Ny],'offsets', [1], 'type_diff','circshift');
CHcirc2 = full(CHcirc'*CHcirc);

%%

miniSS = (maxSS - minSS)*rand(Nx, Ny) + minSS;
miniSS(1,1) = maxSS;
miniSS(end,end) = minSS;

%%
k = [];
scales = 10.^(-2:3);
for ii = 1:length(scales)
        scale = scales(ii);
        T = scale*CH2 + eye(Nx*Ny);%diag(miniSS(:));% + eye(Nx*Ny);
        k(ii) = cond(T);
end

% want to have ratio of I:CH2 to be 10^-0.483 = 0.3289
figure; plot(scales, k)
%%
[V, D] = eig(T);
figure; plot(diag(D))
