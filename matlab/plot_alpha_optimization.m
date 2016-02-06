% plot alpha optimization

clear all;
close all;

set(0, 'DefaultAxesFontSize', 16);

[i r a b1 b2 g] = textread( '../../../experiments/MCM_Annealing_Alpha_Static_Median/parameterFile.txt', '%f %f %f %f %f %f' );

% normalize
nr = 1/max(r)*r;
na = 1/max(a)*a;

f = figure; hold on
plot(i, nr, '-g', i, na, '-r', 'LineWidth', 1)
title('Optimization of alpha', 'fontsize', 16, 'fontweight', 'b')
xlabel('Iteration Step', 'fontsize', 16)
ylabel('Normalized Value', 'fontsize', 16)
h = legend('RMSD FEM', 'alpha', 1);
set(h,'Interpreter','none')
set(h,'FontSize',16);

saveas(f,'../../../images/AlphaOtimizationMCM.eps','epsc2');