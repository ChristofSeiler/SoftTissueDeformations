% draw graphs for documentation
clear all;
close all;

x = [0:300];
y = 1/30*x;

figure; hold on

xlabel('x')
ylabel('spatial displacement')
title('Continuas displacement for inital guess')

plot(x,y)