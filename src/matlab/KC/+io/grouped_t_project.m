function grouped_t_project(filepaths)
%MEMMAP_MOVIE Summary of this function goes here
%   Detailed explanation goes here

if nargin ==1
    file = cell(1,numel(filepaths));
    for i = 1:numel(filepaths)
        [path, name, ext] = fileparts(filepaths(i));
        file{i} = [name ext];
    end
else
    [file, path] = uigetfile({'*.tif'; '*.tiff'},...
                            'Select tif stacks to memory-map',...
                            "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy",...
                            'Multiselect', 'on');                        
end

if ischar(file)
   Y = imread_big(fullfile(path, file)); 
else
    Y = [];             
    for i = 1:numel(file)
        %y = imread_big(fullfile(path, file{i}));
        y = loadtiff(fullfile(path, file{i}));
        Y = cat(3, Y, y);
    end    
end



checked = false;
while checked==false
    slices = input('How many slices in this movie? >> ');
    try
        Y = reshape(Y, size(Y,1), size(Y,2), slices, []);
        checked = true;
    catch
        warning('Invalid- # of frames not divisible by # of slices given.');
    end
end

[d1,d2,d3,T]=matsplit(size(Y));

checked=false;
while checked==false
    tsub = input('Temporal downsampling factor(tsub)? >> ');
    assert(rem(T,tsub)==0, 'User-suggested tsub not valid.');
    checked=true;
end

Yds = zeros(d1,d2,d3,T/tsub,'like',Y);
for z = 1:d3
   y = squeeze(Y(:,:,z,:));
   y = reshape(y,d1,d2,tsub,[]);
   Yds(:,:,z,:) = sum(y,3);
end


[filepath, name, ext] = fileparts(fullfile(path, file));
fname = sprintf('%s_sum%d%s', name, tsub, ext);
fname = fullfile(filepath, fname);

tiffoptions.overwrite = true;
saveastiff(reshape(Yds,d1,d2,[]), fname, tiffoptions);

end
