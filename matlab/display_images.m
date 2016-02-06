clear all;
close all;

I = imread('../DVFRegularization/tissue_image.png');
%L = imread('../DVFRegularization/label_image.png');
G = imread('../DVFRegularization/deformed_grid_image.png');
D = imread('../DVFRegularization/deformed_tissue_image.png');

figure

% Use a 256-value grayscale color map
colormap(gray(256));
ax(1) = subplot(1,3,1);
image(I); title('Undeformed Image')

%ax(2) = subplot(1,3,2);
%image(L); title('Label Image')

ax(2) = subplot(1,3,2);
image(G); title('Deformed Grid Image')

ax(3) = subplot(1,3,3);
image(D); title('Deformed Image')

linkaxes(ax,'xy')
axis(ax,'image')
