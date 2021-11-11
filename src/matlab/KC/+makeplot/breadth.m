function [g, breadth, tuning, unique_groups] = breadth(bin_response, grouping, labels, title_str)
%BREADTH Summary of this function goes here
%   Detailed explanation goes here

[out,idx] = sort(grouping, 'ascend');
bin_response_sorted = bin_response(:,idx);
grouping_sorted = grouping(idx);
labels_sorted = labels(idx);

unique_groups = unique(grouping_sorted);
unique_labels = unique(labels_sorted, 'stable');


fcn = @(x) sum(bin_response_sorted(:, grouping_sorted==x),2);
tuning = cell2mat(arrayfun(fcn, unique_groups, 'UniformOutput', false));

group_counts = arrayfun(@(x) sum(grouping_sorted == x), unique_groups);

cat = 1:max(group_counts);
breadth = cell2mat(arrayfun(@(x) sum(tuning>=x,2), cat, 'UniformOutput', false));


clear g
rcats = repmat(cat, size(breadth,1), 1);
g = gramm('x', breadth(:), 'color', rcats(:));
g.facet_grid(rcats(:), [], 'scale', 'fixed', 'space', 'fixed', 'row_labels', true);
g.stat_bin('geom', 'overlaid_bar', 'edges', -.5:1:10, 'normalization','probability');

g.set_names('x', 'odor breadth', 'row', 'trials');
g.set_text_options('legend_title_scaling', 1);
g.no_legend();

if exist('title_str', 'var') && ~isempty(title_str)
    g.set_title(title_str, 'FontWeight', 'normal');
else
    g.set_title('response breadth', 'FontWeight', 'normal', 'FontSize', 12);
end


g.draw();

end

