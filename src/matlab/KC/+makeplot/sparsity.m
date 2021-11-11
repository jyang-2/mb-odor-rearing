function g = sparsity(response, ti, title_str)
% SPARSITY plots bar graph of odor response sparsities (grouped by odor)
% shows repetitions within one movie
%   Detailed explanation goes here

if ~exist('title_str', 'var') || isempty(title_str)
    title_str = 'sparsity, grouped by odor';
end
K = size(response.bin_response,1);
sparsity = sum(response.bin_response)/K;


%fig = figure('Position', [450 667 350 330]);
codors = ti.pin_odors(ti.si);
%reps = repelem(1:,1,ti.num_stim/3);
uhi = unique(ti.hi);
reps = zeros(size(ti.hi));
for i = 1:numel(uhi)
    idx = (ti.hi==uhi(i));
    reps(idx) = 1:sum(idx);
end

g(1,1) = gramm('x', codors, 'y', sparsity, 'color', reps);
g(1,1).geom_bar('dodge', .6, 'width', .3, 'stacked', false, 'EdgeColor', 'auto');
g(1,1).no_legend();
g(1,1).set_names('x', 'odors', 'y', 'sparsity');
g(1,1).set_color_options('chroma', 0, 'lightness', 0);
g(1,1).axe_property('XTickLabelRotation', 90, 'FontSize', 8, 'Color', 'none');
g(1,1).set_title(title_str, 'FontSize', 10, 'FontWeight', 'normal');
%g.draw();

end

