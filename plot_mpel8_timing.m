truncate = 1;

if truncate
	Nx = 252;
	Ny = 140;
	niter = 2000;
else
	Nx = 256;
	Ny = 144;
	niter = 2000;
end
colors = 'cmbrgk';
xform = @(x) x;
figure; hold on; 
%load(sprintf('./reviv/mpel8_timing_%dx%d_%diter_tunedmu.mat', Nx, Ny, niter), 'time_tri', 'time_tri_max', 'time_alp2c', 'time_tri_2max', 'nrmsd_tri', 'nrmsd_tri_max', 'nrmsd_alp2c', 'nrmsd_tri_2max')
load(sprintf('./reviv/mpel8_timing_%dx%d_%diter_tunedmu.mat', Nx, Ny, niter));
legend_str = {};
plot(cumsum(time_tri), nrmsd_tri,'k'); hold on; legend_str = [legend_str; 'tridiag,40core'];
plot(cumsum(time_tri_max), nrmsd_tri_max,'b'); hold on; legend_str = [legend_str; 'tridiag,80core'];
plot(cumsum(time_tri_2max), nrmsd_tri_2max,'r'); hold on; legend_str = [legend_str; 'tridiag,160core'];
%plot(cumsum(time_tri_mu), nrmsd_tri_mu,'k.'); hold on; legend_str = [legend_str; 'tridiag,40core,tunedmu'];
%plot(cumsum(time_tri_max_mu), nrmsd_tri_max_mu,'b.'); hold on; legend_str = [legend_str; 'tridiag,80core,tunedmu'];
plot(cumsum(time_alp2c), nrmsd_alp2c, 'g'); legend_str = [legend_str; 'AL-P2 circ (to circ x^*)'];

keyboard
for ii = 1:3
if 1
	load(sprintf('./reviv/mpel8_timing_alp2_%dx%d_%dii_tunedmu.mat', Nx, Ny, ii),'xhat_alp2t', 'nrmsd_alp2t', 'time_alp2t');
	if truncate
		plot(cumsum(time_alp2t), nrmsd_alp2t, [colors(ii) '--']); hold on;
	else
		plot(cumsum([time_alp2t(1); time_alp2t]), nrmsd_alp2t, [colors(ii) '--']); hold on;
	end	
else
	plot(cumsum(time_alp2t(:,ii)), nrmsd_alp2t(:,ii), [colors(ii) '--']); hold on;

end

	legend_str = [legend_str; sprintf('AL-P2,%d',ii)];
end
legend(legend_str)
xlabel('wall time (s)')
title(sprintf('NRMSD to x^* over time for [%d %d]', Nx, Ny))
ylabel('NRMSD to x^*');
axis tight
axis([0 110 0 0.14])
axis([1 400 0 0.18]);

%%
figure; hold on; 
%load(sprintf('./reviv/mpel8_timing_%dx%d_%diter_tunedmu.mat', Nx, Ny, niter), 'time_tri', 'time_tri_max', 'time_alp2c', 'time_tri_2max', 'nrmsd_tri', 'nrmsd_tri_max', 'nrmsd_alp2c', 'nrmsd_tri_2max')
load(sprintf('./reviv/mpel8_timing_%dx%d_%diter_tunedmu.mat', Nx, Ny, niter));
legend_str = {};
plot(nrmsd_tri,'k'); hold on; legend_str = [legend_str; 'tridiag,40core'];
plot(nrmsd_tri_max,'b'); hold on; legend_str = [legend_str; 'tridiag,80core'];
plot(nrmsd_tri_2max,'r'); hold on; legend_str = [legend_str; 'tridiag,160core'];
%plot(cumsum(time_tri_mu), nrmsd_tri_mu,'k.'); hold on; legend_str = [legend_str; 'tridiag,40core,tunedmu'];
%plot(cumsum(time_tri_max_mu), nrmsd_tri_max_mu,'b.'); hold on; legend_str = [legend_str; 'tridiag,80core,tunedmu'];
plot(nrmsd_alp2c, 'g'); legend_str = [legend_str; 'AL-P2 circ (to circ x^*)'];

keyboard
for ii = 1:3
if 1
	load(sprintf('./reviv/mpel8_timing_alp2_%dx%d_%dii_tunedmu.mat', Nx, Ny, ii),'xhat_alp2t', 'nrmsd_alp2t', 'time_alp2t');
	if truncate
		plot(nrmsd_alp2t, [colors(ii) '--']); hold on;
	else
		plot(nrmsd_alp2t, [colors(ii) '--']); hold on;
	end	
else
	plot(nrmsd_alp2t(:,ii), [colors(ii) '--']); hold on;

end

	legend_str = [legend_str; sprintf('AL-P2,%d',ii)];
end
legend(legend_str)
xlabel('iteration number')
title(sprintf('NRMSD to x^* over iteration for [%d %d]', Nx, Ny))
ylabel('NRMSD to x^*');
axis tight
axis([0 110 0 0.14])


%set(gca,'DataAspectRatio',[1.35 1 1]);
