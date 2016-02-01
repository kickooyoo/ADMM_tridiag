% experimental data for tridiag_ADMM

tridiag_exp_setup;
niters = 2000;

% load xinfs
%load([home 'Documents/mai_code/ADMM_tridiag/reviv/tri_chcv_5000iter.mat'],'xhat_tri')
load(sprintf('./reviv/x_tri_inf_slice%d_beta%.*d.mat', slice, 3, beta), 'x_tri_inf');
if truncate
%        xhat_tri = reshape(xhat_tri, 256, 144);
        x_tri_inf = x_tri_inf(3:end-2, 3:end-2);
%else
%        x_tri_inf = xhat_tri;
end

% maybe not xinf since not guaranteed to converge
% load([home 'Documents/mai_code/ADMM_tridiag/reviv/al_p2_chcv_5000iter.mat'],'xhat_alp2');
% if truncate
%         xhat_alp2 = reshape(xhat_alp2, 256, 144);
%         x_alp2_inf = xhat_alp2(3:end-2, 3:end-2);
% else
%         x_alp2_inf = xhat_alp2;
% end

%load([home 'Documents/mai_code/ADMM_tridiag/reviv/al_p2_circ_5000iter.mat'],'xhat_alp2');
load(sprintf('./reviv/x_alp2c_inf_slice%d_beta%.*d.mat', slice, 3, beta),'x_alp2c_inf');
if truncate
%        xhat_alp2c = reshape(xhat_alp2, 256, 144);
        x_alp2c_inf = x_alp2c_inf(3:end-2, 3:end-2);
%else
%        x_alp2c_inf = xhat_alp2;
end

if 0
[xhat_alp2c, ~, nrmsd_alp2c, costOrig_alp2c, time_alp2c] = AL_P2_gen(...
        y_noise, F, S, RcircW, SoS, niters, beta, x_alp2c_inf, 'mu', mu, ...
        'zmethod', 'fft', 'inner_iter', 3);

%[xhat_alp2t, ~, nrmsd_alp2t, costOrig_alp2t,time_alp2t] = AL_P2_gen(...
%        y_noise, F, S, RW, SoS, niters, beta, x_tri_inf, 'mu', mu, 'inner_iter', 3);
[xhat_tri, ~, nrmsd_tri, costOrig_tri, time_tri] = tridiag_ADMM(...
        y_noise, F, S, CH, CV, alph, beta, SoS, x_tri_inf, niters,'mu', mu);
save(sprintf('./reviv/mpel8_timing_%dx%d_%diter.mat', Nx, Ny, niters));
end

for ii = 1:10
        [xhat_alp2t, ~, nrmsd_alp2t, costOrig_alp2t, time_alp2t] = AL_P2_gen(y_noise, F, S, R, SoS, niters, beta, x_tri_inf, 'mu', mu,'inner_iter', ii);
        save(sprintf('./reviv/mpel8_timing_alp2_%dx%d_%dii.mat', Nx, Ny, ii),'xhat_alp2t', 'nrmsd_alp2t', 'time_alp2t', 'ii');
end

return;

        [xhat_alp2c, ~, nrmsd_alp2c, costOrig_alp2c, time_alp2c] = AL_P2_gen(...
                y_noise, F, S, Rcirc, SoS, niters, beta, x_alp2c_inf, 'mu', mu, 'zmethod','fft');
[xhat_tri, xsaved_tri, nrmsd_tri, costOrig_tri, time_tri] = tridiag_ADMM(y,F,S,CH,CV,alph,beta,SoS,x_tri_inf,niters,'mu', mu);
return;

[xhat_alp2c, xsaved_alp2c, nrmsd_alp2c, costOrig_alp2c,time_alp2c] = AL_P2_gen(y,F,S,Cdiffs([nx ny],'type_diff','circshift'),SoS,niters,beta,x_alp2c_inf,'mu',mu,'zmethod','fft');
[xhat_alp2, xsaved_alp2, nrmsd_alp2, costOrig_alp2,time_alp2] = AL_P2_gen(y,F,S,R,SoS,niters,beta,x_alp2_inf,'mu', mu);
[xhat_alp2t, xsaved_alp2t, nrmsd_alp2t, costOrig_alp2t,time_alp2t] = AL_P2_gen(y,F,S,R,SoS,niters,beta,x_tri_inf,'mu',mu);

save(sprintf('./reviv/mpel8_timing_%dx%d.mat', Nx, Ny));
send_mai_text('done with mpel8');

return

colors = 'cmbrgk';
figure; hold on; 
legend_str = {};
for ii = 1:5%6
        load(sprintf('./reviv/mpel8_timing_alp2_%dx%d_%dii.mat', 256, 144, ii),'xhat_alp2t', 'nrmsd_alp2t', 'time_alp2t', 'ii');
	plot(cumsum(time_alp2t), nrmsd_alp2t(2:end), [colors(ii) '--']);	
	legend_str = [legend_str; sprintf('AL-P2,%d',ii)];
end
load('./reviv/mpel8_timing_256x144_2000iter.mat')
plot(cumsum(time_tri), nrmsd_tri,'k'); hold on; legend_str = [legend_str; 'tridiag,40core'];
plot(cumsum(time_alp2c), nrmsd_alp2c(2:end), 'g'); legend_str = [legend_str; 'AL-P2 circ (to circ x^*)'];
legend(legend_str)
xlabel('wall time (s)')
title(sprintf('NRMSD to x^* over time for [%d %d]', Nx, Ny))
ylabel('NRMSD to x^*');
axis tight
axis([0 400 0 0.18]);

%set(gca,'DataAspectRatio',[1.35 1 1]);
