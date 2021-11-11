% load from python spks and check out results

S = readNPY('G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\processed_data\2021-03-13\1\Untitled_002\suite2p\plane0\spks.npy');
iscell = readNPY('G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\processed_data\2021-03-13\1\Untitled_002\suite2p\plane0\iscell.npy');
S = S(logical(iscell(:,1)),:);

%% full algorithm
[isort1, isort2, Sm] = mapTmap(S);
imagesc(Sm,[0,3])

%% run map in neurons without smoothing across time sorting
[iclustup, isort, Vout] = activityMap(S);
%%
imagesc(zscore(S(isort,:), 1, 2), [0 3])
