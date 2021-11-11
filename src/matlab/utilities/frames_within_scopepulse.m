function good_idx = frames_within_scopepulse(rT, fT, scope_ict, scope_fct)
%FRAMES_OUTSIDE_SCOPE Returns timepoints within scope pulse
%   [good_frame_times, good_idx] = frames_outside_scope(ti.frameouts.rT, ti.frameouts.fT, ti.scope_ict, ti.scope_fct)
    
    rT = torow(rT);
    fT = torow(fT);
    scope_ict = torow(scope_ict);
    scope_fct = torow(scope_fct);

    assert(numel(rT) == numel(fT));
    assert(numel(scope_ict) == numel(scope_fct));
    
    n_blocks = numel(scope_ict);
    
    frame_len = mean(fT - rT);
    
    good_idx = zeros(1, numel(rT));
    for i = 1:n_blocks
        good_rising = (rT > scope_ict(i)) & (rT < scope_fct(i));
        good_falling = (fT > scope_ict(i)) & (fT < (scope_fct(i) + 0.98*frame_len));
        good_idx(good_rising & good_falling) = 1;
        %idx = find(good_rising & good_falling);
        %good_idx = [good_idx, idx];
    end
    good_idx = logical(good_idx);

end

