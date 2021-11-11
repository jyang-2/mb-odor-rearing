%% Example: SyncData analysis
% Shows how to parse Episode001.h5 file in thorSync folder for frame times
% and stimulus info

close('all'), clear, clc

%% Section 1: set dataset, file info
%h5filepath = fullfile(remyDropboxPath, "mb_odor_rearing/data/2p/2021-04-05/3/SyncData_001)");
%SAVE_DIR = fullfile(remyDropboxPath, "mb_odor_rearing/processed_data/2021-04-05/3/movie_001";
SAVE_DIR = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/processed_data/2021-02-09/1/movie_002";
h5filepath = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/data/2p/2021-02-09/1/SyncData_movie002/Episode001.h5";
thorimagename = 'movie_002';
thorsyncname = 'SyncData_movie002';
date = '2021-02-09';
fly = 1;


%% Section 2 : convert h5 file to downsampled .mat file
if exist(fullfile(SAVE_DIR, 'syncdata.mat'), 'file') 
    % if h5 file already converted
    ai = load(fullfile(SAVE_DIR, 'syncdata.mat')); 
else
    % load and convert file
    ai = thorSync.load_sync_episode(h5filepath);
    k = 30; % downsample by a factor of 30
    ai = thorSync.downsample_ai(ai, k);
    save(fullfile(SAVE_DIR, 'syncdata.mat'), '-struct',  'ai');
end   
%% plot thorSync lines 
figure, 
plot(ai.time, ai.scopePin); hold on; 
plot(ai.time, ai.olfDispPin);
plot(ai.time, ai.Frame_Out);
xlabel('time (sec)');
ylabel('V');
ylim([-1 6])
title('thorsync lines', 'Interpreter', 'none');
legend({'scopePin', 'olfDispPin', 'Frame_Out'}, 'Location', 'best', 'Interpreter', 'none');
legend('boxoff');
%%

exportgraphics(gcf, fullfile(SAVE_DIR, 'syncdata_lines.png'));
%% Section 3: get frame times from thorsync lines
listing = dir('*.tif');
Y = loadtiff(listing.name);
expected_frames = size(Y,3);
fprintf('expected frames : %d\n' , expected_frames);
%%
ai = load('syncdata.mat');
ti = thorSync.ai2ti(ai);

%%
%[frame_times, indexes] = thorSync.generate_frame_time(ai.Frame_Out, ai.time);
[frame_times, indexes] = thorSync.generate_frame_time(ai.FrameOut, ai.time);
ti = struct();
ti.frame_times = frame_times;

%[rLT, fLT, fT] = thorSync.get_frameout_times(ai.Frame_Out', ai.time');
[rLT, fLT, fT] = thorSync.get_frameout_times(ai.FrameOut', ai.time');
ti.frameouts.rT = rLT;   % rise time
ti.frameouts.fT = fLT;   % fall time of frame outs
ti.frameouts.T = (rLT + fLT)/2;


% get timing information about, and number of, blocks and stim
% acquisition block info
[W,ict,fct,~] = pulsewidth(ai.scopePin, ai.time);
ti.num_blocks = length(W);                 % # of blocks
ti.block_len = W;
ti.block_ict = ict;
ti.block_fct = fct;

[W,ict,fct,~] = pulsewidth(ai.olfDispPin, ai.time);  
ti.num_stim = length(W);                  % # of odor pulses
ti.stim_len = W;
ti.stim_ict = ict;
ti.stim_fct = fct;

ti.timepoints = length(ti.frameouts.rT);

% make sure all fields are rows
ti = make_fields_rows(ti);
ti = ti_frame_corrector(ti);
ti.spt = ti.num_stim / ti.num_blocks;               % stimuli per trial

%%

%save(fullfile('ti.mat'), '-struct', 'ti');
save(fullfile('ti0.mat'), '-struct', 'ti');
%save(fullfile(SAVE_DIR, 'ti0.mat'), '-struct', 'ti0');
%%
% ti.baseline_start = arrayfun(@(x) find(T.baseline_idx==x,1,'first'), 1:27);
% ti.baseline_end = arrayfun(@(x) find(T.baseline_idx==x,1,'last'), 1:27);
% ti.peakresp_start = arrayfun(@(x) find(T.peakresp_idx==x,1,'first'), 1:27);
% ti.peakresp_end = arrayfun(@(x) find(T.peakresp_idx==x,1,'last'), 1:27);
% 
% fn = fieldnames(mystruct);
% for k=1:numel(fn)
%     
% end
% 
% ti.block_start = arrayfun(@(x) find(T.block_idx==x,1,'first'), 1:9);
%ti.block_end = arrayfun(@(x) find(T.block_idx==x,1,'last'), 1:9);

%%
frame_times = ti.frame_times;
frame_num = (1:numel(frame_times));

opts = struct('baseline_win', 10, 'peakresp_win', 5, 'trial_win', 20);
%baseline_win = 10;
%peakresp_win = 5;
%trial_win = 20;

olf_pulse_idx = zeros(size(frame_times));
baseline_idx = zeros(size(frame_times));
peakresp_idx = zeros(size(frame_times));
trialtrace_idx = zeros(size(frame_times));

for i = 1:ti.num_stim
   rT = ti.stim_ict(i);
   fT = ti.stim_fct(i);
   olf_pulse_idx(frame_times>=rT & frame_times<=fT) = i;
   baseline_idx(frame_times>=(rT+baseline_win) & frame_times<rT) = i;
   peakresp_idx(frame_times>=rT & frame_times<=rT+peakresp_win) = i;
   trialtrace_idx(frame_times>=(rT+baseline_win) & frame_times<=rT+trial_win)=i;
end

%scope_pulse_idx = repelem((1:ti.num_blocks), 1, (ti.timepoints/ti.num_blocks))
scope_pulse_idx = zeros(size(frame_times));
for i = 1:ti.num_blocks
    rT = ti.block_ict(i);
    fT = ti.block_fct(i);
   scope_pulse_idx(frame_times>=rT & frame_times<=fT+peakresp_win) = i;
end

T = table(frame_num', frame_times', scope_pulse_idx', olf_pulse_idx', baseline_idx', peakresp_idx', trialtrace_idx',...
    'VariableNames', {'frame_num', 'frame_times', 'scope_pulse_idx', 'olf_pulse_idx', 'baseline_idx', 'peakresp_idx', 'trialtrace_idx'});
%%
if all(T.scope_pulse_idx ~= 0)
    fprintf('\tSaving frametimeinfo.csv\n')
    writetable(T, 'frametimeinfo0.csv');
end


