function  memmap_movie(filepaths)
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

checked = false;
while checked==false
    goodz = input('Which slices are good? >> ');
    try
        Y = Y(:,:,goodz,:);
        checked = true;
    catch
        warning('User-suggested slices not valid.');
    end
end

folders = split(path, filesep);
folders = folders(~cellfun(@isempty, folders));
pnam = folders{end};
[savefile, savepath]= uiputfile('*.mat', 'Save memmapped file', fullfile(path, [pnam '.mat']));

Y = single(Y);
sizY = size(Y);
Yr = reshape(Y,prod(sizY(1:end-1)),[]);
nY = min(Yr(:));

savefast(fullfile(savepath, savefile), 'Yr','Y','nY','sizY');

end
