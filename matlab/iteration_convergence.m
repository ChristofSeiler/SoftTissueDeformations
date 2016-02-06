clear all;
close all;

set(0, 'DefaultAxesFontSize', 16);

[ic tc dc ec] = textread( ['time_difference_energy.txt'], '%f %f %f %f' );
[s i] = textread( ['scale_iterations_random.txt'], '%f %f' );

% normalize
ntc = 1/max(tc)*tc;
ndc = 1/max(dc)*dc;
nec = 1/max(ec)*ec;

%scrsz = get(0,'ScreenSize');
%figure('Position', [50 50 scrsz(3)/1.5 scrsz(4)/2.5]); hold on

%subplot(1,2,1);
f1 = figure; hold on
set(gca,'XTick',0:25:250)
plot(ic, ndc, '-g', ic, nec, '--r', 'LineWidth', 3)
title('Convergence', 'fontsize', 16, 'fontweight', 'b')
xlabel('Iteration Step', 'fontsize', 16)
ylabel('Normalized Value', 'fontsize', 16)
h = legend('RMSD Error(FEM,HMRF)', 'Mean Energy Value', '1');
set(h,'Interpreter','none')
set(h,'FontSize',16);
%grid on;
saveas(f1,'../../../images/TwoTissueExpansion/IterationConversion.eps','epsc2');

%subplot(1,2,2);
f2 = figure; hold on
set(gca,'XTick',0:10:100)
plot(ic, tc, '-g', 'LineWidth', 3)
title('Computation Time', 'fontsize', 16, 'fontweight', 'b')
xlabel('Iteration Step', 'fontsize', 16)
ylabel('Time [sec.]', 'fontsize', 16)
%grid on;
saveas(f2,'../../../images/TwoTissueExpansion/ComputationTime.eps','epsc2');

%subplot(1,2,2);
f3 = figure; hold on
set(gca,'XTick',0:1:9)
plot(s, i, '-g', 'LineWidth', 3)
title('Iterations Per Scale Level', 'fontsize', 16, 'fontweight', 'b')
xlabel('Scale Level', 'fontsize', 16)
ylabel('Iterations', 'fontsize', 16)
%grid on;
saveas(f3,'../../../images/TwoTissueExpansion/ScaleIterations.eps','epsc2');
