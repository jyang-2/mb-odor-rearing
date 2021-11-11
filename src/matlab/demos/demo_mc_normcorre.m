PRJ_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing";

get_mov_dir = @(x,y,z) fullfile(PRJ_DIR, x, char(y), z);
get_tiff_dir = @(x,y,z,blk) fullfile(PRJ_DIR, x, char(y), z, char(blk));

%-------------------------------------
%              parameters
%-------------------------------------
flydate = '2018-08-30';
flynum = 2;
mov = 'movie_001';
block = [0, 1, 2];

params = {};

%fly 1
params{end+1} = {'2018-08-30', '1', 'movie_002', 0};
params{end+1} = {'2018-08-30', '1', 'movie_002', 1};
params{end+1} = {'2018-08-30', '1', 'movie_002', 2};

params{end+1} = {'2018-08-30', '1', 'movie_003', 0};
params{end+1} = {'2018-08-30', '1', 'movie_003', 1};
params{end+1} = {'2018-08-30', '1', 'movie_003', 2};

% fly 2
params{end+1} = {'2018-08-30', '2', 'movie_001', 0};
params{end+1} = {'2018-08-30', '2', 'movie_001', 1};
params{end+1} = {'2018-08-30', '2', 'movie_001', 2};

params{end+1} = {'2018-08-30', '2', 'movie_002', ''};

% fly 3
params{end+1} = {'2018-08-30', '3', 'movie_001', 0};
params{end+1} = {'2018-08-30', '3', 'movie_001', 1};
params{end+1} = {'2018-08-30', '3', 'movie_001', 2};

params{end+1} = {'2018-08-30', '3', 'movie_002', 0};
params{end+1} = {'2018-08-30', '3', 'movie_002', 1};
params{end+1} = {'2018-08-30', '3', 'movie_002', 2};

% fly 4
params{end+1} = {'2018-08-30', '4', 'movie_001', []};

params{end+1} = {'2018-08-30', '4', 'movie_002', []};

params{end+1} = {'2018-08-30', '4', 'movie_003', []};

for i = 1:numel(params)
   prm = params(i,:); 
   MOV_DIR = char(get_tiff_dir(prm{1:3}));
   if prm{end} == []
       TIFF_DIR = MOV_DIR;
   else
      TIFF_DIR = char(get_tiff_dir(prm));
   end

end


% load xml_meta.json
xml_meta = struct();

% load tiff files if not saved as .mat
load_tiff_files = ~exist(fullfile(TIF_DIR, 'stk.mat'), 'file');

if load_tiff_files
    tiff_files = glob(fullfile(TIFF_DIR, 'stk_*.tif'));
    Y = [];
    for i = 1:numel(tiff_files)
       file = tiff_files{i};
       disp(file);
       y = loadtiff(file); 
       [d1, d2, frames] = size(y);   
       Y = cat(3, Y, y);
    end    
    % Y = reshape(Y, d1, d2, d3, []);
    savefast(fullfile(TIFF_DIR, 'stks.mat'), 'Y');
    clear Y y tiff_files;
end

% get dimensions from matfile
fname_mmap = fullfile(TIFF_DIR, 'stks.mat');
mm = matfile(fname_mmap);
[d1,d2,d3,T] = size(mm, 'Y');


%%
close all, clear, clc
TIFF_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/data/processed_data/2021-08-30/2/movie_001/1";
TIFF_DIR = char(TIFF_DIR);
%%
tiff_files = glob(fullfile(TIFF_DIR, 'stk_*.tif'));
z = 16;

Y = [];
for i = 1:numel(tiff_files)
   file = tiff_files{i};
   disp(file);
   y = loadtiff(file); 
   [d1, d2, frames] = size(y);   
   Y = cat(3, Y, y);
end
Y = reshape(Y, d1, d2, z, []);
T = size(Y, 4);
siz = size(Y);

savefast(fullfile(TIFF_DIR, 'stks.mat'), 'Y');

clear Y y;
%%
fname_mmap = fullfile(TIFF_DIR, 'stks.mat');
mm = matfile(fname_mmap);
[d1,d2,d3,T] = size(mm, 'Y');

% close parallel pool
delete(gcp('nocreate'))

% run moco
options_nr = NoRMCorreSetParms('d1', d1,...
    'd2', d2,...
    'd3',d3,...
    'grid_size', [64, 64, 3],...
    'overlap_pre', [8,8,1],...
    'overlap_post', [32 32 4],...
    'mot_uf',[4,4,1],...
    'bin_width', 50,...
    'max_shift',[32 32 2],...
    'max_dev',5,...
    'init_batch', 200,...
    'us_fac',4,...
    'iter', 2,...
    'shifts_method', 'FFT',...
    'upd_template', true,...
    'h5_filename', fullfile(TIFF_DIR, 'motion_corrected_nr.h5'),...
    'output_type', 'hdf5'...
    );
[M_nr,shifts_nr,template_nr] = normcorre_batch_even(Y, options_nr);
%%



%% save to tiff stacks

chunk_size = 323;

fid = 1;
fprintf(fid, '==================================================================\n');
fprintf(fid, 'SAVING MOTION CORRECTED MOVIES\n');
fprintf(fid, '...................................................................\n');
fprintf(fid, 'tiff directory = %s\n', TIFF_DIR);
fprintf(fid, '...................................................................\n');

fprintf(fid, 'original size = %s\n', strjoin(string(size(M2)), ", "));
for t0 = 1:chunk_size:T
    t1 = min(t0+chunk_size-1, T);
    m = M2(:,:,:,t0:t1);
    savename = sprintf('moco_nr_%04d.tif', floor(t0/chunk_size));
    fprintf(fid, '\t%s ==> size: [%s]', savename, strjoin(string(size(m)), ", "));
    
    m = reshape(m, d1, d2, []);    
    saveastiff(m, char(fullfile(TIFF_DIR, savename)), struct('overwrite', true));    
    fprintf(fid, '.......done saving\n')
end
%%

drop_planes = [1, 2, 16];




img = imtile(imadjustn(mat2gray(template2)));
saveastiff(single(img), fullfile(TIFF_DIR, 'tiled_template.tif'), struct('overwrite', true));
saveastiff(single(template2), fullfile(TIFF_DIR, 'template_nr.tif'), struct('overwrite', true));
%%
%[cY,mY,vY] = motion_metrics(fname_mmap, [], [], 'Y');
for t0 = 1 %:chunk_size:T
    t1 = min(t0+chunk_size-1, T);
    m = M2(:,:,:,t0:t1);
    [cM2,mM2,vM2] = motion_metrics(m);
    
end



%%

save(fullfile(TIFF_DIR, 'moco_nr.mat'), 'options_nr', 'shifts2', 'template2'); % 'cM2', 'mM2', 'vM2', 'cY', 'vY', );