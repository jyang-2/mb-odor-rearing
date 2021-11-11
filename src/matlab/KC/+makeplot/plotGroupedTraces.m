function fig = plotGroupedTraces(trial_traces,cids, grouping,labels, stim_frame,title_str)
%PLOTGROUPEDTRACES Summary of this function goes here
% INPUTS
%  trial_traces : [time x num_stim x K]
%      grouping : [1 x num_stim]
%        labels : [1 x num_stim], string array or cellstr
%    stim_frame : frame # within trial when stimulus turns on
%      title_str: title
%
% OUTPUTS
%            fig: figure with traces plotted


%% 
yspacing = 5;
xspacing = 10;
scaling = 1;

%%
[stim_len, num_stim, K] = size(trial_traces);

[grouping_sorted,idx] = sort(grouping, 'ascend');
labels_sorted = labels(idx);

unique_groups = unique(grouping_sorted);
unique_labels = unique(labels_sorted, 'stable');
num_groups = numel(unique_groups);
%%
%hi = unique(ti.hi);
%cids = (1:30)+0;

dy = 1:numel(cids);
dy = (dy-1) * -yspacing;

dx = 1:num_groups;
dx = (dx-1) * (stim_len + xspacing);

%frame_offset = (0:(num_stim-1))*xspacing;

ts = 1:stim_len;

%%

fig = figure('renderer', 'painter');
hold on


% plot traces
for i = 1:numel(unique_groups)
    for j = 1:numel(cids)
        traces = trial_traces(:,grouping==unique_groups(i),cids(j));

        lh = plot(ts+dx(i), traces+dy(j), 'k-');
        set(lh, 'Color', [0 0 0 .7]);
        
        text(-25, dy(j), num2str(cids(j)), 'FontSize', 8, 'HorizontalAlignment', 'right');
        
    end
    set(gca, 'ColorOrderIndex', 1);
end

% draw stimulus on lines
stim_on = stim_frame + dx;
llstim = arrayfun(@(x) xline(x,':'), stim_on);
xticks(stim_on)

% odor labels
xticklabels(unique_labels);
xtickangle(90);
set(gca, 'XAxisLocation', 'top', 'FontSize', 8);

% title
if exist('title_str', 'var') && ~isempty(title_str)
    title(title_str, 'FontWeight', 'normal', 'FontSize', 10, 'Interpreter', 'none');
end

end

