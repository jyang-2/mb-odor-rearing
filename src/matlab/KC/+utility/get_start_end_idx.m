function [vals, start_idx, end_idx] = get_start_end_idx(vec, vals)
%GET_START_END_IDX Summary of this function goes here
%   Detailed explanation goes here

    arguments
        vec (1,:) {isnumeric} 
        vals = unique(vec(vec~=0 & ~isnan(vec)))
    end
    
    start_idx = arrayfun(@(x) find(vec==x,1,'first'), vals);
    end_idx = arrayfun(@(x) find(vec==x,1,'last'), vals);
    
end

