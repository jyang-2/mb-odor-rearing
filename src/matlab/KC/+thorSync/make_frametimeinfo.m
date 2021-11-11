function T = make_frametimeinfo(frame_times, ti, opts)
%MAKE_FRAMETIMEINFO Summary of this function goes here
%   Detailed explanation goes here

frame_num = (1:numel(frame_times));

baseline_win = opts.baseline_win;
peakresp_win = opts.peakresp_win;
trial_win = opts.trial_win;

% baseline_win = -10;
% peakresp_win = 5;
% trial_win = 20;

scope_pulse_idx = thorSync.get_frametime_idx(frame_times, ti.block_ict, ti.block_fct);
olf_pulse_idx = thorSync.get_frametime_idx(frame_times, ti.stim_ict, ti.stim_fct);

baseline_idx = thorSync.get_frametime_idx(frame_times, ti.stim_ict+baseline_win, ti.stim_ict);
peakresp_idx = thorSync.get_frametime_idx(frame_times, ti.stim_ict, ti.stim_ict + peakresp_win);
trialtrace_idx = thorSync.get_frametime_idx(frame_times, ti.stim_ict+baseline_win, ti.stim_ict+trial_win);

%baseline_idx = zeros(size(frame_times));
%peakresp_idx = zeros(size(frame_times));
%trialtrace_idx = zeros(size(frame_times));


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


end

