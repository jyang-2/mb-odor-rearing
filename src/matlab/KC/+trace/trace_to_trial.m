function [C_df_trial] = trace_to_trial(C_df, ti)
%GET_TRIAL_TRACES Returns time series split into trials
%                 dim : [# cells x time x trial]
%
%  C_df_trial = trace_to_trial(C_df, ti)

[K,T] = size(C_df);
assert(rem(T,ti.num_trials)==0,...
    'Number of frames cannot be evenly split into desired # of trials.');
try
    trace.get_interval_vec(ti);
catch
    error('Timing of events and/or trial length inconsistent across trials.');
end

C_df_trial = reshape(C_df', [], ti.num_trials, K);  % time x trials x cells
C_df_trial = permute(C_df_trial, [3 1 2]);    % cells x time x trial

end



