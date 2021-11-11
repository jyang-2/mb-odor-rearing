function [fig, hax] = plot_trial_factors(M, R, trials, c)
import utility.*

%PLOT_TRIAL_FACTORS plot_trial_factors(M, R, trial, trial_labels)
%   plots trial factors, in stimulus order and grouped by trial type
%   c: colors
neuron_factors = M.U{1};
time_factors = M.U{2};
trial_factors = M.U{3};
lambda = M.lambda;

[ustim, ~, trial_color_idx] = unique(trials);
cmap_keys = viridis(numel(ustim));
cmap_values = viridis(numel(ustim));
cmap_dict = containers.Map(keySet,valueSet,'UniformValues', isUniform);

% sort trials
[sorted_trials,idx] = sort(trials);
num_stim = numel(trials);
rep = utility.get_rep(trials);

if ~exist('trial_labels','var')
    trial_labels = string(1,num_stim);
end


%hax = gobjects(R,2);
hax = gobjects(R,4);
rloc = utility.line_loc(rep);
lloc = utility.line_loc(sorted_trials);


fig = figure('Visible','off');
tl = tiledlayout(R,4, 'padding', 'compact', 'tilespacing', 'compact');


for i = 1:R    
    % left column, 
    hax(i,1) = nexttile;
    nfac = neuron_factors(:,i);
    bar(nfac, 'FaceColor', 'k', 'EdgeColor', 'k');
    U = M.U{1};
    xl = [0 size(M,1)+1];
    yl = [min( 0, min(U(:)) ), max( 0, max(U(:)) )];
    xlim(xl);
    ylim(yl);
    ylabel(num2str(i));

    
    hax(i,2) = nexttile;
    tsfac = time_factors(:,i);
    plot(tsfac, '-k');
    U = M.U{2};
    xl = [0 size(M,2)+1];
    yl = [min( 0, min(U(:)) ), max( 0, max(U(:)) )];
    xlim(xl);
    ylim(yl);
    
    % trials in stimulus order
    hax(i,3) = nexttile;
    tfac = trial_factors(:,i);
    scatter(hax(i,3), 1:num_stim, tfac, trial_color_idx);
    %vline(rloc,'k:');
   
    % trials sorted by odor
    hax(i,4) = nexttile;
    scatter(hax(i,4), 1:num_stim, tfac(idx)', trial_color_idx(sort_idx));
    
    %vline(lloc,'k:');
    %arrayfun(@(x) xline(x), lloc);
    
end

trialax = hax(:,3:4);

% set(trialax,'YLim', [0 1]);
set(hax(1:end-1,:), 'XTickLabel',[], 'Color', 'none', 'Box', 'on');
set(hax(R,:), 'Color', 'none', 'Box', 'on');
% set(trialax(1,2), 'XTick', utility.label_loc(trials), 'XAxisLocation', 'top', 'FontSize', 8);
% xticklabels(trialax(1,2), unique(trial_labels(idx),'stable'));
% xtickangle(trialax(1,2), 90);

hax(1,1).Title.String = 'neurons';
hax(1,2).Title.String = 'time';
hax(1,3).Title.String = 'trials';
hax(1,4).Title.String = 'trials (sorted)';

xline([ 0; lloc], '-', unique(trial_labels(idx),'stable'));

xx = cell2mat(get(trialax, 'XLim'));
for i = 1:R
    text(trialax(i,1), xx(i,2), 1, ['\lambda=' num2str(lambda(i),'%.2f') ' '],...
        'FontSize', 8, 'Interpreter', 'tex',...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    %line(trialax(i,1), [rloc; rloc] ,[0 1],'Color', 'k', 'LineStyle', ':');
    %line(trialax(i,2), [lloc; lloc] ,[0 1],'Color', 'k', 'LineStyle', ':');
    
end

tt = sprintf('Trial factors (%d)', R);
tl.Title.String = tt;
tl.Title.FontSize = 10;
fig.Visible = 'on';




end

