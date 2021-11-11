function frametimeinfo = load_frametimeinfo(folder)
%LOAD_FRAMETIMEINFO Loads table of timing information from
%                   frametimeinfo.csv
% Inputs
%   - folder
%
% Outputs
%   - frametimeinfo: table w/ variables {'frame_num', 'frame_times',
%                      'scope_pulse_idx', 'olf_pulse_idx', 'baseline_idx', 
%                       'peakresp_idx','trialtrace_idx'}
frametimeinfo = readtable(fullfile(folder, 'frametimeinfo.csv'));
end

