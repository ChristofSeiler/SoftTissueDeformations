% plot fem convergence

clear all;
close all;

[e t d] = textread( 'tumor_elements_time_displacement.txt', '%f %f %f' );

% normalize
nt = 1/max(t)*t;
nd = 1/max(d)*d;

f = figure; hold on
%plot(e, nt, '-g', e, nd, '-r', 'LineWidth', 1)
plot(e, d, '-g', 'LineWidth', 1)
title('FEM convergence', 'fontsize', 16, 'fontweight', 'b')
xlabel('Iteration step', 'fontsize', 16)
ylabel('Normalized value', 'fontsize', 16)
h = legend('time', 'displacement', 1);
set(h,'Interpreter','none')
set(h,'FontSize',16);

saveas(f,'../../../images/AlphaOtimizationMCM.eps','epsc2');