clear all;
close all;

set(0, 'DefaultAxesFontSize', 16);

[pva xva yva] = textread( '../DVFRegularization/v_plot_abaqus.txt', '%f %f %f' );
[pv xv yv] = textread( '../DVFRegularization/v_plot_9.txt', '%f %f %f' );
[pha xha yha] = textread( '../DVFRegularization/h_plot_abaqus.txt', '%f %f %f' );
[ph xh yh] = textread( '../DVFRegularization/h_plot_9.txt', '%f %f %f' );
[pd1a xd1a yd1a] = textread( '../DVFRegularization/d1_plot_abaqus.txt', '%f %f %f' );
[pd1 xd1 yd1] = textread( '../DVFRegularization/d1_plot_9.txt', '%f %f %f' );
[pd2a xd2a yd2a] = textread( '../DVFRegularization/d2_plot_abaqus.txt', '%f %f %f' );
[pd2 xd2 yd2] = textread( '../DVFRegularization/d2_plot_9.txt', '%f %f %f' );

f = figure; hold on
%set(gca,'YTick',-50:25:50)
plot(pva, yva, '-g', pv, yv, '--r', 'LineWidth', 3)
title('Vertical Profile', 'fontsize', 16, 'fontweight', 'b')
xlabel('Pixel', 'fontsize', 16)
ylabel('Displacement [pixel]', 'fontsize', 16)
h = legend('ABAQUS', 'HMRF', 2);
set(h,'Interpreter','none')
set(h,'FontSize',16);
saveas(f,'../../../images/TwoTissueExpansion/ProfileVertical.eps','epsc2');

f = figure; hold on
%set(gca,'YTick',-50:25:50)
plot(pha, xha, '-g', ph, xh, '--r', 'LineWidth', 3)
title('Horizontal Profile', 'fontsize', 16, 'fontweight', 'b')
xlabel('Pixel', 'fontsize', 16)
ylabel('Displacement [pixel]', 'fontsize', 16)
h = legend('ABAQUS', 'HMRF', 1);
set(h,'Interpreter','none')
set(h,'FontSize',16);
saveas(f,'../../../images/TwoTissueExpansion/ProfileHorizontal.eps','epsc2');

for i=1:size(pd1a)
   ma1(i) = norm([xd1a(i) yd1a(i)]');
   m1(i) = norm([xd1(i) yd1(i)]');
end
f = figure; hold on
%set(gca,'YTick',-50:25:50)
plot(pd1a, ma1, '-g', pd1, m1, '--r', 'LineWidth', 3)
title('Diagonal Profile', 'fontsize', 16, 'fontweight', 'b')
xlabel('Pixel', 'fontsize', 16)
ylabel('Displacement [pixel]', 'fontsize', 16)
h = legend('ABAQUS', 'HMRF', 1);
set(h,'Interpreter','none')
set(h,'FontSize',16);
saveas(f,'../../../images/TwoTissueExpansion/ProfileDiagonal1.eps','epsc2');

for i=1:size(pd2a)
   ma2(i) = norm([xd2a(i) yd2a(i)]');
   m2(i) = norm([xd2(i) yd2(i)]');
end
f = figure; hold on
%set(gca,'YTick',-50:25:50)
plot(pd2a, ma2, '-g', pd2, m2, '--r', 'LineWidth', 3)
title('Diagonal Profile', 'fontsize', 16, 'fontweight', 'b')
xlabel('Pixel', 'fontsize', 16)
ylabel('Displacement [pixel]', 'fontsize', 16)
h = legend('ABAQUS', 'HMRF', 2);
set(h,'Interpreter','none')
set(h,'FontSize',16);
saveas(f,'../../../images/TwoTissueExpansion/ProfileDiagonal2.eps','epsc2');
