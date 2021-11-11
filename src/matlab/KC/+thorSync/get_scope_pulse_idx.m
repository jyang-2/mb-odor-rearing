function scope_pulse_idx = get_scope_pulse_idx(frame_times, scope_ict, scope_fct)
%GET_SCOPEP Summary of this function goes here
%   Detailed explanation goes here


nBlocks = numel(scope_ict);

scope_pulse_idx = zeros(size(frame_times));
for i = 1:nBlocks
    rT = scope_ict(i);
    fT = scope_fct(i);
   scope_pulse_idx(frame_times>=rT & frame_times<=fT) = i;
end

end

