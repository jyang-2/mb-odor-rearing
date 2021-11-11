function fname = downsample_movie(memmapped_filepath,tsub)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
[~,~,ext] = fileparts(memmapped_filepath);
if strcmp(ext, '.mat')
    load(memmapped_filepath, 'Y');    
elseif strcmp(ext, '.tif') || strcmp(ext, '.tiff')
    tinf = tiffinfo(memmapped_filepath);
    d1 = tinf.Height;
    d2 = tinf.Width;
    d3 = tinf.ParsedImageDescription.slices;
    T = tinf.ParsedImageDescription.frames;
    sizY = [d1 d2 d3 T];
    
    Y = imread_big(memmapped_filepath);
    Y = reshape(Y,sizY);
end

Y = downsample_data(Y,'time',tsub,1);
Y = uint16(Y);
sizY = size(Y);
Yr = reshape(Y,prod(sizY(1:end-1)),[]);
nY = min(Yr(:));

[filepath, name, ext] = fileparts(memmapped_filepath);
fname = sprintf('%s_tsub%d.mat', name, tsub);
fname = fullfile(filepath, fname);

savefast(fname, 'Yr','Y','nY','sizY');
end

