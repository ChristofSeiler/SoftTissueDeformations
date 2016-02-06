clear all;
close all;

set(0, 'DefaultAxesFontSize', 16);

% read data
[x1 x2 e] = textread( '../EnergyFunctionPlot/x_y_energy_value.txt', '%f %f %f' );

% figure; hold on
% title('Energy Field with 2 x 1 x 1 image')
% xlabel('pixel 1')
% ylabel('pixel 2')
% zlabel('energy')
% plot3(x1, x2, e, 'o')

% fill matrix
for i = 1:size(e)
    M(x1(i)+31,x2(i)+31) = e(i);
end

pos = [50, 50, 672, 504];

f = surf(M)
title('Energy Field on Two Pixel Image', 'fontsize', 16, 'fontweight', 'b')
xlabel('Pixel 1')
ylabel('Pixel 2')
zlabel('Energy')
colorbar
view(128,27)

saveas(f,'../../../images/EnergyFunctionPlot.eps','epsc2');