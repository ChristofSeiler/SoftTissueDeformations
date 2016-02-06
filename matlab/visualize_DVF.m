clear all; 
close all;

% [x f] = textread( 'DisplacementMagnitude.rpt', '%f %f' );
% figure; hold on
% plot(x, f)
% 
% x = [0:19];
% [f] = textread( 'Energy.rpt', '%f' );
% figure; hold on
% plot(x, f)

[x f] = textread( '../vtkL2Interpolate/energyFile.txt', '%f %f' );
figure; hold on
plot(x, f)
