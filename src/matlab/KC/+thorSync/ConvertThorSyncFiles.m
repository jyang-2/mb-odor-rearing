
function [] = ConvertThorSyncFiles(folder, varargin)

default_outfolder = fullfile(folder, "thorSync");

p = inputParser;
validFolderName = @(x) ischar(x) || isStringScalar(x);
addRequired(p,'folder',@isfolder);
addOptional(p,'outfolder', default_outfolder, validFolderName);
parse(p,folder,varargin{:});

folder = p.Results.folder;
outfolder = p.Results.outfolder;
   
listing = dir(fullfile(folder, "*/Episode001.h5"));
disp(folder)

if ~exist(outfolder, 'dir')
    mkdir(outfolder)
end

%%
for i = 1:numel(listing)
   filepath = fullfile(listing(i).folder, listing(i).name);
   
   S = thorSync.FnLoadSyncEpisode(filepath);
   
   folders = regexp(filepath, filesep, 'split');
   thorSyncName = folders{end-1};
   fname = strcat(thorSyncName, ".mat");
   str = sprintf('%s saved', fullfile(outfolder, fname));
   save(fullfile(outfolder, fname), '-struct', 'S');
end


