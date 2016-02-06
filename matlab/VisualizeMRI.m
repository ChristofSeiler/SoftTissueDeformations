load mri
fid = fopen('126_1_20.img','r');
out = fread(fid,[256,256],'int16');
%readanalyze('IBSR_01_ana');
image(out)
axis image
colormap(map)

