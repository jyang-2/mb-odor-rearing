function ti = frametimeinfo2ti(frametimeinfo, trial_list)
%FRAMETIMEINFO2TI Takes a frametimeinfo table and computes trial frame
%information for calculating response vectors
%   
% inputs:
%   - frametimeinfo : table w/ variables {'frame_num', 'frame_times',
%                      'scope_pulse_idx', 'olf_pulse_idx', 'baseline_idx', 
%                       'peakresp_idx','trialtrace_idx'}
%   - trial_list: list of which stimuli (represented in 'olf_pulse_idx')
%                     to include
%                  if empty, the function just uses all the stimuli
%
% outputs 
%    - ti: named for old trial information struct, but here it just contains 
%             frame/indexing info for calculating response vectors           
%       - baseline_start
%       - baseline_end
%       - peakresp_start
%       - peakresp_end

    
arguments
   frametimeinfo table
   trial_list = unique(frametimeinfo.olf_pulse_idx(frametimeinfo.olf_pulse_idx > 0))'
end

disp(trial_list)
% build struct containing frame timing info
ti = struct();
ti.baseline_start = torow(arrayfun(@(x) find(frametimeinfo.baseline_idx==x, 1, 'first'), trial_list));
ti.baseline_end = torow(arrayfun(@(x) find(frametimeinfo.baseline_idx==x, 1, 'last'), trial_list));
ti.peakresp_start = torow(arrayfun(@(x) find(frametimeinfo.peakresp_idx==x, 1, 'first'), trial_list));
ti.peakresp_end = torow(arrayfun(@(x) find(frametimeinfo.peakresp_idx==x, 1, 'last'), trial_list));

end

