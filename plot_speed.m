load tri_chcv_5000iter
load tri_wavelet_chcv_5000iter
load al_p2_rcirc_500iter
load al_p2_wavelet_rcirc_500iter
%
figure; plot(cumsum([0; time_alp2]), nrmsd_alp2);
hold on;
mask = generate_mask('slice67',1,nx,ny);   
nrmsd1_tri = calc_NRMSE_over_mask(SoS(:),x_tri_inf(:),mask);  
nrmsd1_triw = calc_NRMSE_over_mask(SoS(:),x_triw_inf(:),mask);
%plot(cumsum([0; time_tri])-(0:niters)'.*0.001053, [nrmsd1_tri nrmsd_tri],'r');
%plot(cumsum([0; time_triw])-(0:niters)'.*0.001053, [nrmsd1_triw nrmsd_triw],'g'); 
%plot(cumsum([0; time_tri])+(0:niters)'.*0.003159, [nrmsd1_tri nrmsd_tri],'r');
%plot(cumsum([0; time_triw])+(0:niters)'.*0.003159, [nrmsd1_triw nrmsd_triw],'g'); 
plot(cumsum([0; time_tri]), [nrmsd1_tri nrmsd_tri],'k--');

%%
close all
figure; plot(cumsum([0; time_alp2w]),nrmsd_alp2w,'b','LineWidth',2); hold on;
plot(cumsum([0; time_triw]), [nrmsd1_triw nrmsd_triw],'r','LineWidth',2);
% trick to add symbol
%plot(sum([0; time_alp2w(1:50)]), nrmsd_alp2w(51),'bs','LineWidth',2);
%plot(sum([0; time_triw(1:100)]), nrmsd_triw(100),'ro','LineWidth',2);
axis([0 50 0 0.2]);
title('Speed Comparison of Algorithms','FontSize',16);
hlegend = legend('AL-P2',['ADMM-tridiag' char(10) '(emulated wall time)']);
set(hlegend,'FontSize',12);
xlabel('time (seconds)','FontSize',16);
ylabel('NRMSD to x^(*^)','FontSize',16);



