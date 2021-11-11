function loc = line_loc(grouping)
%LINE_LOC Summary of this function goes here
%   Detailed explanation goes here

[grouping_sorted, sort_idx] = sort(grouping, 'ascend');

if iscategorical(grouping)
    [~, ~, ustim_num] = unique(grouping);
    grouping_sorted = ustim_num(sort_idx);
end
% find locations to draw grid lines for grouping odors
loc = find(diff(grouping_sorted))+.5;

end

