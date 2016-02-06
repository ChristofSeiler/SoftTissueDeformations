clear all;
close all;

% read data
[p9 x9 y9] = textread( '../DVFRegularization/h_plot_9.txt', '%f %f %f' );

xAxis='horizontal cross section [pixel]';
yAxis='displacement [pixel]';

% comparison with analitic solution
pos = [50, 50, 672, 504];
figure('Position',pos); hold on
title('level 9')
xlabel(xAxis)
ylabel(yAxis)
%plot(p9,-x9,'-rs')
% analytic solution for 3 tissue
% YM: 80, 20, 10
%y10 = piecewise_3_tissues(x9);
%x10_1=0:1:101; x10_2=100:1:201; x10_3=201:1:302;
%plot(p9, -x9,'-r', x10_1, x10_1.*2.75/101, '--b', x10_2, 2.75+x10_1.*5.5/101, '--m', x10_3, 8.25+x10_1.*22/101, '--g', 'LineWidth', 2)
%h = legend('MRF approximation', 'Young''s modulus = 80 MPa (analytic solution)', 'Young''s modulus = 20 MPa (analytic solution)', 'Young''s modulus = 10 MPa (analytic solution)', 2);
% YM: 40, 20, 10
%x10_1=0:1:101; x10_2=100:1:201; x10_3=201:1:302;
%plot(p9, -x9,'-r', x10_1, x10_1.*4.3/101, '--b', x10_2, 4.3+x10_1.*8.6/101, '--m', x10_3, 12.9+x10_1.*17.2/101, '--g', 'LineWidth', 2)
%h = legend('MRF approximation', 'Young''s modulus = 40 MPa (analytic solution)', 'Young''s modulus = 20 MPa (analytic solution)', 'Young''s modulus = 10 MPa (analytic solution)', 2);
% analytic solution for 2 tissue
% YM: 80, 10
%x10_1=0:1:201; x10_2=201:1:302; x10_3=0:1:101;
%plot(p9, -x9,'-r', x10_1, x10_1.*3/101, '--b', x10_2, 6+x10_3.*24/101, '--g', 'LineWidth', 2)
%h = legend('MRF approximation', 'Young''s modulus = 80 MPa (analytic solution)', 'Young''s modulus = 10 MPa (analytic solution)', 2);
% YM: 20, 10
x10_1=0:1:150; x10_2=151:1:301;
plot(p9, -x9,'-r', x10_1, x10_1.*10/151, '--b', x10_2, 10+x10_1.*20/151, '--g', 'LineWidth', 2)
h = legend('MRF approximation', 'Young''s modulus = 20 MPa (analytic solution)', 'Young''s modulus = 10 MPa (analytic solution)', 2);
set(h,'Interpreter','none')
grid on;