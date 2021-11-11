clc; clear variables; close all;
config = ReadYaml('config.yaml');
meta = ReadYaml('metadata.yaml');
PRJ_DIR = config.lab.prjdir;
BASE_DIR = fullfile(PRJ_DIR, 'processed_data', '2021-03-13', '1');
%%
rng(42);

S = readNPY("/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/processed_data/2021-03-13/1/Untitled_002/suite2p/plane0/spks.npy");
iscell = readNPY("/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/processed_data/2021-03-13/1/Untitled_002/suite2p/plane0/iscell.npy");
ops = readNPY("/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/processed_data/2021-03-13/1/Untitled_002/suite2p/plane0/ops.npy");
%%
S = S(iscell(:,1) == 1, :); % only use cells

% zscore
for i = 1:size(S,1)
    S(i,:) = (S(i,:) - mean(S(i,:)))./std(S(i,:));
end

S = max(-4, min(8, S)) + 4; % Bound zscores on -4 to 8
S = S./12; % Bound values on 0 to 1
%%
efilt = exp(- linspace(0,50,200) / (ops.tau * ops.fs));
%efilt = efilt/sum(efilt);
%sout = nan(389,3829);
sout = nan(size(spks));
for i = 1:height(spks)
    y = spks(i,:);
    sout(i,:) = conv(y, efilt,'same');
end
sout = zscore(sout,1,2);
%%
S = sout(iscell(:,1) == 1, :); % only use cells

% embed from 2904 dimensions to 1 dimension
S_embed = tsne(S, 'NumDimensions', 1, 'Algorithm', 'exact', ...
    'Standardize', false,'Perplexity', 50, 'Verbose', 2);

% sort embedding
[~, index] = sort(S_embed, 'descend');
S = S(index, :);

% smoothing of data (optional)
for i = 1:size(S,1)
    S(i,:) = smooth(S(i,:));
end

% plot
threshold = 0.3; % do not show cell activity < threshold

%cmap = colormap('bone');
%cmap = cmap(end:-1:1,:);
%colormap(cmap);
figure()
imagesc(S); axis xy; box off;
caxis([threshold 1]);
xlabel('Time (seconds)');
ylabel('Cells');
set(gca,'YTick',0:10:10000);