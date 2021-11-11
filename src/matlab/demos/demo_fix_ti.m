%SAVE_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/processed_data/2021-03-23/1/movie";
SAVE_DIR = "G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\processed_data\2021-03-23\1\movie";

frametimeinfo = io.load_frametimeinfo(SAVE_DIR);
ti = io.load_ti(SAVE_DIR);
single_expt_metadata = io.load_single_expt_metadata(SAVE_DIR);
olfactometer = single_expt_metadata.olfactometer;
stimulus = single_expt_metadata.stimulus;

chanA = getChannelInfo(stimulus.channelA, olfactometer);
chanB = getChannelInfo(stimulus.channelB, olfactometer);

stim_labels = strings(25,1);
tbl = table(chanA.hi', chanB.hi', chanA.odor', chanA.conc', chanB.odor', chanB.conc',...
    'VariableNames', {'hi1', 'hi2', 'odorA', 'concA', 'odorB', 'concB'});
for i = 1:height(tbl)
    stim_labels(i) = sprintf('%s [%d], %s [%d]', tbl.odorA(i), tbl.concA(i), tbl.odorB(i), tbl.concB(i));
end
tbl.stim_labels = stim_labels;


tbl_stimtype = readtable("G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\processed_data\2021-03-23\1\movie\stimtype_lookup.xlsx");
tbl = join(tbl, tbl_stimtype);
disp(tbl)

stim = struct();
stim.stim_idx = 1:27;
stim.block_idx = repelem((1:ti.num_blocks),1, ti.spt);
[~, stim.stim_on, stim.stim_off] = utility.get_start_end_idx(frametimeinfo.olf_pulse_idx);
[~, stim.block_start, stim.block_end] = utility.get_start_end_idx(frametimeinfo.scope_pulse_idx, stim.block_idx);
[~, stim.peakresp_start, stim.peakresp_end] = utility.get_start_end_idx(frametimeinfo.peakresp_idx);
[~, stim.baseline_start, stim.baseline_end] = utility.get_start_end_idx(frametimeinfo.baseline_idx);
%%
ti.stim = stim;
save(fullfile(SAVE_DIR,'ti.mat'),'-struct','ti');

