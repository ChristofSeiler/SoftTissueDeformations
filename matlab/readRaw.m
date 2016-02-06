function  rawData = readRaw(filename,imSize,type,endian,skip);

% Reads a raw file
% Inputs: filename, image size, image type, byte ordering, skip
% If you are using an offset to access another slice of the volume image
% be sure to multiply your skip value by the number of bytes of the 
% type (ie. float is 4 bytes).
% Inputs: filename, image size, pixel type, endian, number of values to
% skip.
% Output: image


fid = fopen(filename,'rb',endian);
if (fid < 0)
    fprintf('Filename %s does not exist\n',filename);
    rawData = -1;
else
    status = fseek(fid,skip,'bof');
	if status == 0
        rawData = fread(fid,prod(imSize),type);
        fclose(fid);
        if (length(imSize == 3))
            slices = length(rawData)/imSize(1)/imSize(2);
            imSize(3) = slices;
        end
        rawData = reshape(rawData,imSize);
	else
        rawData = status;
	end
end