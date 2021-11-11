function h = imMat(response_mat, grouping, labels, title_str)
%CORRMAT - plots heatmap of response_mat, 
%   Detailed explanation goes here


[out,idx] = sort(grouping, 'ascend');
response_mat_sorted = response_mat(:,idx);
response_mat_sorted = response_mat_sorted(idx,:);
grouping_sorted = grouping(idx);
labels_sorted = labels(idx);

unique_groups = unique(grouping_sorted);

% find locations to draw axis labels
fcn = @(x) mean(find(grouping_sorted==x));
label_loc = arrayfun(fcn, unique_groups);

% get labels for each group
fcn = @(x) find(grouping_sorted==x,1);
group_labels = labels_sorted(arrayfun(fcn, unique_groups));

% find locations to draw grid lines for grouping odors
line_loc = find(diff(grouping_sorted))+.5;

% calculate correlation coefficients
%R = corrcoef(response_mat_sorted, 'Rows', 'pairwise');             % presentation order

% get and set colormap
load('favorite_colormaps.mat', 'BlueRed');
colormap(BlueRed);

% draw heatmap and colorbar
h=imagesc(response_mat_sorted, [-1 1]);
h.AlphaData = ones(size(h.CData)); 
h.AlphaData(isnan(h.CData)) = 0;
axis image; 
cb = colorbar;
    
% draw gridlines
vline(line_loc, 'k-');
hline(line_loc, 'k-');

% adjust appearence of labels, text
set(gca, 'XTick', label_loc, 'XTickLabel', group_labels, 'XAxisLocation', 'top',...
       'YTick', label_loc, 'YTickLabel', group_labels, 'TickLabelInterpreter', 'none');
xtickangle(90); 
ax = gca;
ax.YAxis.FontSize=8;
ax.XAxis.FontSize=8;
ax.TickDir = 'out';

% set colorbar font size
cb.FontSize = 10;

% display title
if exist('title_str', 'var') && ~isempty(title_str)
    title(title_str, 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'none');
end

end
