colors = 'cmbrgkym';
start_ndx = 2;
end_ndx = niters + 1;
xform = @(x) 20*log10(x);
figure; hold on;
legend_str = {};
for ii = 1:length(nthread_vals);
	plot(cumsum(time_tri(start_ndx:end_ndx,ii)), nrmsd_tri(start_ndx:end_ndx,ii), [colors(ii) '--+']); 
	legend_str = [legend_str; sprintf('%d threads', nthread_vals(ii))];
end
legend(legend_str)
xlabel('wall time (s)')
title(sprintf('NRMSD to x^* over time for [%d %d]', Nx, Ny))
ylabel('NRMSD to x^* (dB)');
axis tight




