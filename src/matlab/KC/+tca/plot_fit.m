function fig = plot_fit(nfac, Err,Sim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

num_comp = size(Err,1);

fig = figure('Renderer', 'painters', 'Units', 'inches');

tl = tiledlayout(1,2);
tl.Padding = 'compact';
tl.TileSpacing = 'compact';

m = size(Err, 2);
A = 0.05;
xx = repmat(nfac',1,m);
xx = xx + rand(size(Err)) * 0.25;

pltvars = {'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'Marker', 'o', 'LineWidth', 0.7};

% scree plot of error
nexttile;
%lh = plot(nfac, Err, pltvars{:});
lh = plot(xx(:), Err(:), 'ko', pltvars{:});
ylabel('model error')
xlabel('# factors');
ylim([0 1]);
title('error (scree plot)')
pbaspect(gca, [1 1 1])
set(gca, 'Color', 'none');

% plot of model similarity
nexttile;
%lh = plot(nfac, Sim, 'ko', pltvars{:});
lh = plot(xx(:), Sim(:), 'ko', pltvars{:});
ylabel('model similarity')
xlabel('# factors');
ylim([0 1]);
title('similarity');
pbaspect(gca, [1 1 1]);
set(gca, 'Color', 'none');

fig.Position(3:4) = [6 3];
end
