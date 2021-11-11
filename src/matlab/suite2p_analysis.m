%file = "G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\data\processed_data\2021-08-30\2\movie_001\0\suite2p\combined\TENSOR_TRIALS.mat";

%PRJ_DIR = "G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rear
%PRJ_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing";

if isunix
    PRJ_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing";
elseif ispc
    PRJ_DIR = "G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing";
end
PROC_DIR = fullfile(PRJ_DIR, 'data', 'processed_data');

file_pfo = fullfile(PRJ_DIR, 'data', 'processed_data', '2021-06-16', '1', '**/plane0', '**/stat.npy');
listing = dir(file_pfo);
%%

statfile = fullfile(listing(2).folder, listing(2).name);
statfolder = fileparts(statfile);

% load data from cell
Fall = load(fullfile(statfolder, 'Fall.mat'));
iscell = Fall.iscell(:,1);
cellprob = Fall.iscell(:,2);

TENSOR_TRIALS = load(fullfile(statfolder, 'TENSOR_TRIALS.mat'));

% fix df_stimulus, make odors categorical
ntmeta = load(fullfile(statfolder, "nt_meta_data.mat"));
ntmeta.meta_yaml = jsondecode(ntmeta.meta_yaml);
%%
zdfof_tensor = TENSOR_TRIALS.zdfof_tensor;
dfof_tensor = TENSOR_TRIALS.dfof_tensor;
trial_ts = TENSOR_TRIALS.trial_ts;

[ranked_cells, cellsort_idx] = sort(cellprob, 'descend');
zdfof_tensor_ranked = zdfof_tensor(cellsort_idx, :, :);
%%
stackedplot(squeeze(zdfof_tensor(1,:,:))');
%%
% fix df_stimulus
df_stimulus = struct2table(ntmeta.df_stimulus);
for field = ["odor_a", "odor_b", "stimname"]
    df_stimulus.(field) = string(df_stimulus.(field));
end
valueset = ["pfo", "1-5ol", "VA", "1-6ol", "EP", "1-6ol+EP"];
df_stimulus.stimnames = categorical(df_stimulus.stimname,valueset,'Ordinal',true);
[sorted_stimnames, stimname_idx] = sort(df_stimulus.stimname);
%%
Fz = zscore(Fc, [], 2);
F_low = lowpass(Fz',0.2,3.0)';
%% 
stim_ict = ntmeta.sync_meta.stim_ict(TENSOR_TRIALS.stim_idx);

stim_times = [];
rep = [];


stim_times = cell(1, numel(stim_ict));
for i = 1:numel(stim_ict)
    ict = stim_ict(i);
    disp(ict)
    
    stim_times{i} = find((ntmeta.timestamps>ict+0.5) & (ntmeta.timestamps<ict+5));
end
%%
svm_table = [];
for i = 1:18
    df = table(stim_times{i}');
    df{:,2} = repmat(df_stimulus.stimname(i), numel(stim_times(i)), 1);
    svm_table = [svm_table; df];
end


%% Linear SVM
stimname = df_stimulus.stimname;
X = Fgood(:, svm_table{:,1}); 
Y = (svm_table{:,2}=='1-6ol') | svm_table{:,2}=="1-6ol+EP";
Lambda = logspace(-6,-0.5,11);
rng(10); % For reproducibility
CVMdl = fitclinear(X,Y,'ObservationsIn','columns','KFold',5,...
    'Learner','logistic','Solver','sparsa','Regularization','lasso',...
    'Lambda',Lambda,'GradientTolerance',1e-8);
%%
Mdl1 = CVMdl.Trained{1};
numCLModels = numel(CVMdl.Trained);
ce = kfoldLoss(CVMdl);

Mdl = fitclinear(X,Y,'ObservationsIn','columns',...
    'Learner','logistic','Solver','sparsa','Regularization','lasso',...
    'Lambda',Lambda,'GradientTolerance',1e-8);
numNZCoeff = sum(Mdl.Beta~=0);
%%
figure;
[h,hL1,hL2] = plotyy(log10(Lambda),log10(ce),...
    log10(Lambda),log10(numNZCoeff)); 
hL1.Marker = 'o';
hL2.Marker = 'o';
ylabel(h(1),'log_{10} classification error')
ylabel(h(2),'log_{10} nonzero-coefficient frequency')
xlabel('log_{10} Lambda')
title('Test-Sample Statistics')
hold off
%%
rng('default') % For reproducibility
X = Fgood(stim_times, :);
Y = ones(size(X));
%%
Fc = Fall.F - 0.7 * Fall.Fneu;
[n_cells, T] = size(Fc);

Fc = cumsum(Fc, 2);
F = diff(X,1,2)
Fc = reshape(Fc', T/blocks, n_cells, blocks);
disp(size(Fc)); % time x cells x blocks

Fc = pagetranspose(Fc);

fcn = @(block_struct) block_struct.data(:,:,[2 1 3]);
M = movsum(A,3,'Endpoints','discard')

%%


ttrials = load(file);
zdfof_tensor = ttrials.zdfof_tensor;
zdfof_tens
disp(size(zdfof_tensor))
%%
zdfof_tensor = permute(zdfof_tensor, [1, 3, 2]);
disp(size(zdfof_tensor));
zdfof_tensor = zdfof_tensor(iscell, :, :);
%%
zdfof_tensor_sorted = zdfof_tensor(:,:,stimname_idx);
show_cells = 1:200;

figure()
tiledlayout(1, 18);
for i = 1:height(df_stimulus)
   nexttile;
   imagesc(zdfof_tensor_sorted(show_cells,:,i), [0, 2]);
   xline(find(ttrials.trial_ts==0), 'w-', sstimname(i), 'LabelOrientation', 'horizontal')

end
%%

model = cp_als(tensor(zdfof_tensor), 10); % fit CP model with 10 components
visualize_neuron_ktensor(model) % produces a nice plot for you
%%
Respvec = load(fullfile(folder, 'RESPVEC.mat'));

show_cells = 1:200;

peak_amp_sorted = Respvec.peak_amp(:, stimname_idx);
peak_amp_sorted = zscore(peak_amp_sorted, [], 1);
figure, imagesc(peak_amp_sorted(show_cells, :), [0 3]);
%%



%%



%%
listing = dir(fullfile(fileparts(folder),  "plane6", "reg_tif", "*.tif"));
I = []
for i = 1:length(listing)
   I = cat(3, I, loadtiff(fullfile(listing(i).folder, listing(i).name)));
end
%%
[d1, d2, T] = size(I);
%%
stim_on_frame = floor(interp1(ntmeta.timestamps, 1:T, df_stimulus.stim_ict));
fs = 3;
peak_start = 1.0;
peak_len = 1.0

baseline_im = cell(1, height(df_stimulus));
peak_im = cell(1, height(df_stimulus));
baseline_len = 10;


for istim = 1:numel(df_stimulus.stim_ict)
    
    baseline_frames = floor(stim_on_frame(istim)-baseline_len*fs):floor(stim_on_frame(istim)-1);
    peak_frames = ceil(stim_on_frame(istim)+1*fs) : floor(stim_on_frame(istim)+(1+5)*fs);
    
    baseline_im{istim} = mean(I(:,:,baseline_frames), 3);
    peak_im{istim} = mean(I(:,:,peak_frames), 3);
end
%%


%%
saveastiff(single(baseline_im), 'baseline.tif', struct('overwrite', true));
saveastiff(single(peak_im), 'peak.tif', struct('overwrite', true));
    

%%
odor_order = ["pfo", "1-5ol", "VA", "1-6ol", "EP", "1-6ol+EP"];
[n_cells, n_trials, n_time] = size(cells_x_trials_x_time);

%mean_trials = cell(1, numel(odor_order));


odor_list = odor_order(ismember(odor_order, df_stimulus.stimname));
n_odors = length(odor_list);
mean_trials = cell(1,n_odors);

tl = tiledlayout(1,n_odors,'TileSpacing','compact');
cids = 1:200

for i = 1:n_odors
   trials = cells_x_trials_x_time(:, df_stimulus.stimname == odor_order(i), :);
   mean_trials{i} = squeeze(nanmean(trials, 2));
   
   nexttile;
   imagesc(squeeze(mean_trials{i}(cids,:)))
    colormap(gray)
    xline(find(trial_ts==0), 'w-', odor_list(i), 'LabelOrientation', 'horizontal')
end
mean_trials = cat(3, mean_trials{:});
%%
montage(mean_trials)
%%

%%
