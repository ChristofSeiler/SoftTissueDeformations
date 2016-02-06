I = imread('../../../images/SphincterAbaqusElasticRegistration.png');
J = imread('../../../images/SphincterElasticRegistration.png');
K = imabsdiff(I,J);
imshow(K,[]) % [] = scale data automatically

I = imread('../../../images/SphincterAbaqusHMRF.png');
J = imread('../../../images/SphincterHMRF.png');
K = imabsdiff(I,J);
%imshow(K,[]) % [] = scale data automatically
