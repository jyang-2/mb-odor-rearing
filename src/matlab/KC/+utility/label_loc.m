function loc = label_loc(grouping)
%LABEL_LOC Summary of this function goes here
%   Detailed explanation goes here

grouping_sorted = sort(grouping, 'ascend');
unique_groups = unique(grouping);

% find locations to draw axis labels
fcn = @(x) mean(find(grouping_sorted==x));
loc = arrayfun(fcn, unique_groups);
end
