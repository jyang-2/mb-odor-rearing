clear;
clc;
gcp;

% same demo as demo_script.m but using the class @CNMF

%% load file
base_nam = "movie_002";
nam = 'fn_001';
date='2021-04-05';
fly=2;

PROC_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/processed_data/2021-04-05/2/movie_002";
tifflist = dir("G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\processed_data\2021-04-05\2\movie_002\cnmf\reg_tif\*.tif");
%%
for i = 2

%     ti.baseline_start = arrayfun(@(x) find(T.baseline_idx==x,1,'first'), 1:27);
%     ti.baseline_end = arrayfun(@(x) find(T.baseline_idx==x,1,'last'), 1:27);
%     ti.peakresp_start = arrayfun(@(x) find(T.peakresp_idx==x,1,'first'), 1:27);
%     ti.peakresp_end = arrayfun(@(x) find(T.peakresp_idx==x,1,'last'), 1:27);
    
    disp(tifflist(i))
    nam = sprintf('fn_%s', num2str(i,'%03d'));
    filename = fullfile(tifflist(i).folder, tifflist(i).name);
    Y = loadtiff(filename);
    [d1,d2,T]=size(Y);
    if ~isfolder(num2str(i))
        mkdir(num2str(i));
    end

    CNM.loadArray(Y);
    CNM.optionsSet(options);                        % setup the options structure
    CNM.gSig = options.gSig;
    %% Process the dataset
    close all
    CNM.preprocess();             % preprocessing (compute some quantities)
    CNM.initComponents(K);      % initialization
    CNM.plotCenters()           % plot center of ROIs detected during initialization
        F = getframe(gca);
        imwrite(F.cdata, fullfile(num2str(i), 'plotInitCenters.png'));

    %% Update components
    CNM.updateSpatial();        % update spatial components
    CNM.updateTemporal(0);      % update temporal components (do not deconvolve at this point)
    %% Plot contours
    close all
    CNM.cm = [];
    CNM.contours=[];
    CNM.COM();
    CNM.K = size(CNM.A,2);
    CNM.CI=mY;
    CNM.correlationImage();

    
    CNM.plotContours(1);
end
% %%
% is_memmaped=false;
% 
% nam = 'fn_002';
% filename = '../movie_002_s2pReg.tif';
% disp('=================================================================')
% fprintf('\nfilename = %s\n', filename);
% 
% mat_filename = sprintf('%s_cnmf2d.mat', nam);
% mat_filepath = fullfile('cnmf', mat_filename);
% fprintf('\nmat_filepath = %s\n', mat_filepath);
% 
% base_mat_filename = sprintf('%s_cnmf.mat', nam);
% base_mat_filepath = fullfile(base_mat_filename);
% %%
% 
% Y = loadtiff(filename,1,730);
% [d1,d2,T]=size(Y);
%%

save('cnmf.mat','base_nam', 'nam', 'date', 'fly','options', 'CNM', 'CNMc','S','ti0');

%%
CNM = CNMF;
K = 400;        % # of expected components
tau = [3 3];        % 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = .85;                                  % merging threshold

options = CNMFSetParms(...   
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'fr', 4.5,...
    'p',p,...                                   % order of AR dynamics    
    'gSig',tau,...                                % half size of neuron
    'merge_thr',merge_thr,...                        % merging threshold  
    'nb', 2,...                                  % number of background components    
    'min_SNR',3,...                             % minimum SNR threshold
    'noise_norm', false,...
    'rem_prct', 0,...
    'bas_nonneg', false,...
    'space_thresh',0.4,...                      % space correlation threshold
    'time_thresh',0.4,...                      % space correlation threshold
    'cnn_thr',0.2,...                            % threshold for CNN classifier    
    'init_method', 'greedy',...
    'search_method', 'ellipse',...
    'dist', 2,...
    'max_size', 4,...
    'min_size', 3, ...
    'clos_op', strel('square',2),...
    'spatial_method', 'regularized',...
    'thr_method', 'nrg',...
    'nrgthr', .98,...
    'min_pixel', 10,...
    'cluster_pixels', true,...
    'extract_max', true,...
    'bd', 30, ...
    'fast_merge', true ...
    );
%'extract_max', true
%'dist', 3,...
%'max_size', 3,...
%'min_size', 1,...

%%
% Below is the standard processing pipeline. This processing can be
% executed in one shot using the CNM.fit function:
 %CNM.fit(Y,options,K)
%% load the dataset and create the object
%CNM.readFile('../fn_001.tif');                         % insert path to file here  
CNM.loadArray(Y);
CNM.optionsSet(options);                        % setup the options structure
CNM.gSig = options.gSig;
%% Process the dataset
close all
CNM.preprocess();             % preprocessing (compute some quantities)
CNM.initComponents(K);      % initialization
CNM.plotCenters()           % plot center of ROIs detected during initialization
%%
sn = reshape(CNM.sn,CNM.dims);
%Ysd = std(CNM.Y,1,3);
CNM.CI = sn;
CNM.manuallyRefineComponents();
%% Update components
CNM.updateSpatial();        % update spatial components
CNM.updateTemporal(0);      % update temporal components (do not deconvolve at this point)
%% Plot contours
close all
CNM.cm = [];
CNM.contours=[];
CNM.COM();
CNM.K = size(CNM.A,2);
mY = mean(CNMc.Y,3);
CNM.CI=mY;
%CNM.correlationImage();

figure()
CNM.plotContours(1);

%% Plot contours
CNMc  = copy(CNM);
CNMc.cm = [];
CNMc.contours=[];
CNMc.COM();
CNMc.K = size(CNM.A,2);
CNMc.CI=[];
%CNM.correlationImage();
figure
CNMc.plotContours();

%% component classification
CNM.evaluateComponents();   % evaluate spatial components based on their correlation with the data
CNM.eventExceptionality();  % evaluate traces

ind_throw = find(~CNM.keep_eval);
ind_keep = find(CNM.keep_eval);
%%
ind_keep = true(size(CNM.A,2),1);
ind_keep(ind_throw) = false;
ind_keep = find(ind_keep);
%%
figure
CNM.contours=[];
CNM.plotContours(1,[],[],ind_throw); % plot components to get rid of
%CNM.plotContours(1);
%%
%get rid of components if sure
CNM.keepComponents(ind_keep);      % keep the components that are above certain thresholds
%% merge components
CNM.merge();                            %% merge found components
CNM.displayMerging();

%% repeat processing and extract dF/F values
CNM.updateSpatial();
CNM.updateTemporal(0);
%% plot components in GUI
CNMc.plotComponentsGUI();     % display all components
%%
figure()
tl = tiledlayout(1,8);
linecount = 0;
for i=1:400
    if mod(linecount,50)==0
        axis tight
        nexttile();
    end
    plot(CNMc.C_df(i,:)-i);
    hold on
    linecount = linecount+1;
    
end
axis tight;
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
%%
responseOptions.std_thr = 3;
responseOptions.thr = 0.5;
responseOptions.max_k = 1;
response = compute_response_vectors(S.sDFF, ti0, responseOptions);
%%
ti0 = struct();
ti0.scope_pulse_start = find(T.scope_pulse_idx==2,1,'first');
dframe = ti0.scope_pulse_start-1;
ti0.num_stim = 3;
ti0.block_idx = ti.block_idx(ti.block_idx==2);
ti0.baseline_start = ti.baseline_start(ti.block_idx==2)-dframe;
ti0.baseline_end = ti.baseline_end(ti.block_idx==2)-dframe;
ti0.peakamp_start = ti.peakresp_start(ti.block_idx==2)-dframe;
ti0.peakamp_end = ti.peakresp_end(ti.block_idx==2)-dframe;
%%
stim = table(repmat("1-5ol",3,1),conc1(ti.block_idx==2)', repmat("1-6ol",3,1),conc2(ti.block_idx==2)', 'VariableNames', {'odor1', 'conc1', 'odor2', 'conc2'})

%%
g = makeplot.sparsity(response, ti0);
figure, g.draw();

set(gca, 'Color', 'none'); 

annotation('textbox', [0.7 0.85 0.15 .1], 'FitBoxToText','on',...
            'String', response.options.str,...
            'FontSize', 8, 'EdgeColor', 'none', 'BackgroundColor', 'none');
%%
fig_folder = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/2019 data analysis/2019-07-19_analysis/2/cnmf/movie002_rig(0)_bg_2D/figures";
F = getframe(gca);
imwrite(F.cdata, fullfile(fig_folder, 'movie001_z01_contours.png'));
%%
v_dims = [224 224 10];
v_d = prod(v_dims);
zvol = zeros([224 224 10]);
zvol(:,:,z) = z;
zvol = logical(reshape(zvol, numel(zvol), 1));

vA = zeros(v_d, size(CNMc.A,2));
vA(zvol,:) = CNMc.A;
%%
vA = sparse(vA);
%%

%%
Atemp = CNMc.A;
for i = 1:K
    Atemp(:,i) = Atemp(:,i) * neuron_sn(i);
end
%spatial = sum(CNMc.A,2);
spatial = sum(Atemp,2);
spatial = reshape(spatial, 224,224);
%%
C_df = CNMc.C_df;
neuron_sn = GetSn(C_df);
% for i = 1:K
%    sn = GetSn 
% end
    
%%
figure
imagesc(spatial*10)
