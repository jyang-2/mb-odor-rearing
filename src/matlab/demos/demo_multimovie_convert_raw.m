
listing = dir('**/Image_001_001.raw');

%listing = listing(~isroimask,:);
%disp({listing.folder}')
folders = {listing.folder}';
rel_folders = cellfun(@(x) erase(x, pwd), folders, 'UniformOutput', false);
%%
PRJ_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing";
SAVE_DIR = "G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\processed_data\2021-08-24\2";

MOVIE_DIR = fullfile(PRJ_DIR, 'data', '2p', '2021-08-30', num2str(2), 'movie_002');
PROC_DIR = fullfile(PRJ_DIR, 'data', 'processed_data',...
                    relativepath(MOVIE_DIR, fullfile(PRJ_DIR, 'data', '2p')))
nFolders = numel(listing);
for i = 5:nFolders
    folder = listing(i).folder;
    [imData,ThorImageExperiment] = thorImage.load_experiment_raw(folder);
    
    [filepath,name,ext] = fileparts(rel_folders{i});
    savename = fullfile(SAVE_DIR, filepath, name, [name '.tif']);
    saveastiff(uint16(imData), char(savename), struct('overwrite', true));
    disp(savename)
end

%%
flydate = '2021-08-30';
flynum = 2;
movnam = 'movie_002';

PRJ_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing";

MOVIE_DIR = fullfile(PRJ_DIR, 'data', '2p', '2021-08-30', num2str(2), movnam);
PROC_DIR = fullfile(PRJ_DIR, 'data', 'processed_data',...
                    relativepath(MOVIE_DIR, fullfile(PRJ_DIR, 'data', '2p')));
                
ti = load(fullfile(PROC_DIR, 'ti.mat'));
siz = [256, 256];
[imData,ThorImageExperiment] = thorImage.load_experiment_raw(MOVIE_DIR);
imData = reshape(imData, siz(1), siz(2), 20, []);
imData = imData(:,:,1:16,:);
%%
imData = reshape(imData, siz(1), siz(2), 16, [], ti.num_blocks);
%%
for iblock = 1:ti.num_blocks
   y  = squeeze(imData(:, :, :, :, iblock));
   y = reshape(y, 256, 256, []);
   
   savename = sprintf('stk_%04d.tif', iblock);
   savename = fullfile(PROC_DIR, savename);
   saveastiff(uint16(y), char(savename), struct('overwrite', true));
   disp(savename)
end
