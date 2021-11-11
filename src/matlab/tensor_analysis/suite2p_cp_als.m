close all, clear, clc

if isunix
    PRJ_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing";
elseif ispc
    PRJ_DIR = "G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing";
end
PROC_DIR = fullfile(PRJ_DIR, 'data', 'processed_data');

stat_listing = dir(fullfile(PROC_DIR, '2021-08-30', '**/combined', '**/stat.npy'));

%%
% load stat.npy file from suite2p folder
statfile = fullfile(stat_listing(1).folder, stat_listing(1).name);
suite2p_folder = fileparts(statfile);

% load original suite2p data
Fall = load(fullfile(suite2p_folder, 'Fall.mat'));

% correct neuropil, extract cell probs
% try to work with full fluorescence matrix
Fc = Fall.F - 0.7 * Fall.Fneu;
iscell = Fall.iscell(:, 1);
cellprob=Fall.iscell(:, 2);
iplane = readNPY(fullfile(suite2p_folder, 'iplane.npy'));

[n_cell, T] = size(Fc);

%%

ipix = sub2ind([yL(iplane) xL(iplane)], yext, xext);

img{iplane}(ipix, 1) = rand;
img{iplane}(ipix, 2) = 1;
img{iplane}(ipix, 3) = 1;
%% fix df_stimulus, make odors categorical

% load metadata
ntmeta = load(fullfile(statfolder, "nt_meta_data.mat"));
ntmeta.meta_yaml = jsondecode(ntmeta.meta_yaml);

% convert to table
df_stimulus = struct2table(ntmeta.df_stimulus);
for field = ["odor_a", "odor_b", "stimname"]
    df_stimulus.(field) = string(df_stimulus.(field));
end

% make df_stimname.stimnames ordinal
valueset = ["pfo", "1-5ol", "VA", "1-6ol", "EP", "1-6ol+EP"];
df_stimulus.stimnames = categorical(df_stimulus.stimname,valueset,'Ordinal',true);
[stimname_sorted, stimname_idx] = sort(df_stimulus.stimname);
%%
cmap_odor_set = turbo(6);
[~, ia] = ismember(stimname_sorted, valueset);
[~, ib] = ismember(df_stimulus.stimname, valueset);

cmap_presentation = cmap_odor_set(ia, :);
cmap_sorted = cmap_odor_set(ia, :);

%% load TENSOR_TRIALS from suite2p results folder

TENSOR_TRIALS = load(fullfile(statfolder, 'TENSOR_TRIALS.mat'));
zdfof_tensor = TENSOR_TRIALS.zdfof_tensor;
zdfof_tensor = permute(zdfof_tensor, [1 3 2]);

dfof_tensor =  TENSOR_TRIALS.dfof_tensor;
dfof_tensor = permute(dfof_tensor, [1 3 2]);

%% Run cp_als model
data = tensor(dfof_tensor);

close all;

R = 10;
n_fits = 1;
err = zeros(n_fits,1);
for n = 1:n_fits
    % fit model
    est_factors = cp_als(tensor(data),R, 'tol', 1e-6);
    
    if n==1
        first_factors = est_factors;
    end
    
    % store error
    err(n) = norm(full(est_factors) - data)/norm(data);
    
    
    % visualize fit for first several fits
    if n > 1
        % score aligns the cp decompositions
        [sc, est_factors] = score(est_factors, first_factors);
        
        % plot the estimated factors
        viz_ktensor(est_factors, ... 
            'Plottype', {'bar', 'line', 'scatter'}, ...
            'Modetitles', {'neurons', 'time', 'trials'})
        set(gcf, 'Name', ['estimated factors - fit #' num2str(n)])
   end
end

%[~, ~, ic] = unique(df_stimulus.stimnames);
%[fig, hax] = tca.plot_trial_factors(est_factors, R, ic, stimname);

info = viz_ktensor(est_factors, ...
'Plottype', {'bar', 'line', 'scatter'}, ...
'Modetitles', {'neurons', 'time', 'trials'},...
'Plotcolors',{'k', 'k', 'k'});
set(info.FactorAxes(1:end-1,:), 'Color', 'none', 'Box', 'on');
set(info.FactorAxes(end,:), 'Color', 'none', 'Box', 'on');
%%
neuron_factors = est_factors.U{1};
time_factors = est_factors.U{2};
trial_factors = est_factors.U{3};
lambda = est_factors.lambda;

%%

Y = tsne(neuron_factors);
scatter(Y(:,1), Y(:,2));
%%
figure
tiledlayout(1,4)

nexttile;
imagesc(trial_factors(sort_idx, :)' ); 

nexttile([1, 3]);
imagesc(time_factors'), colormap('turbo');
%%
figure();
imagesc(neuron_factors(cids, :));
linkdata on

%% Process traces
Fc_split = mat2cell(Fc,n_cell, (T/3)*ones(1,3));

Fc_lowpass_block = cell(1, 3);
qt_block = cell(1, 3);
Fc_lowpass_qt_block = cell(1, 3);

for i = 1:3
    Fc_lowpass_block{i} = lowpass(Fc_split{i}', 0.2, 3)';
    qt_block{i} = quantile(Fc_lowpass_block{i}, 0.6, 2) - quantile(Fc_lowpass_block{i}, 0.05, 2);
    Fc_lowpass_qt_block{i} = Fc_lowpass_block{i}./qt_block{i};
    %qt_block{i} = quantile(Fc_split{i}, 0.7, 2);
    %Fc_qt_block{i} = Fc_split{i}./qt_block{i};
end
%%
