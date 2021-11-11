dirman = DirectoryManager('datestr', '2018-10-21', 'fly', 1);
dirman.findCodeDir;
dirman.findAnalysisDir;
dirman.findSaveDir;
dirman.savedir = fileparts(dirman.savedir);
%%
hi = ti.hi;
tuning = respan.get_tuning(response, hi);
breadth = sum(tuning.fractrials > 0, 2);
nodors = sum(tuning.fractrials>0,2);

Y = tsne(response.peakAmp);
%%
fig_tsne = tsneChart('xData', Y(:,1), 'YData', Y(:,2), 'group', nodors);
%%
fig_tsne.TitleText = [db.tag ': tsne (color coded by resp. breadth'];
%%
fig_tsne.drawRoi()
%%
figs = fig_tsne.plotRoiTraces(S.DFF, S.sDFF, ti, response);
%%
utility.print_figures([fig_tsne.Parent figs], 'folder', dirman.savedir,...
    'base_plot_name', '002_tsneBreadth',...
    'idstr', dirman.datefly);
%%

addRequired(p, 'figs')
addParameter(p, 'folder', tempdir);
addParameter(p, 'base_plot_name', 'plots');
addParameter(p, 'idstr', '');
addParameter(p, 'vstr', datetime('now', 'Format', fmt));
addParameter(p, 'filefmts', {'-dpdf', '-dsvg'});
%%
figure(1)
clf;
gs = gscatter(Y(:,1), Y(:,2), breadth)
ax = gca;

title('tsne of peakAmp');
ax.Legend.Title.String = 'breadth';
ax.Legend.Title.FontWeight = 'normal';
ax.Legend.Color = 'none';

%%
fhs = findobj(ax, 'Type', 'images.roi.freehand');
for i=1:numel(fhs)
   fhs(i).Label = num2str(i); 
end
%%
figs = gobjects(0);

traces_per_fig = 100;
for i = 1:numel(fhs)
    roi = fhs(i);
    pts_selected = inROI(roi, Y(:,1), Y(:,2));
    cids = find(pts_selected);
    
    while(numel(cids)>100)
        ft = CNMFPP.plotOverlaidTraces(S.DFF, S.sDFF, cids(1:traces_per_fig), ti, response);
        title(i);    
        cids(1:traces_per_fig)=[];
        figs(numel(figs)+1) = ft;
        
    end
    ft = CNMFPP.plotOverlaidTraces(S.DFF, S.sDFF, cids, ti, response);
    title(i);
    figs(numel(figs)+1) = ft;
end
%%
dar = get(gca, 'DataAspectRatio')

hax = findobj(figs, 'Type', 'axes')
set(hax, 'DataAspectRatio', dar);
%%
parentdir = fullfile(analysis_2019, 'Multifly analysis')

%%
savedir = char("/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/2019 data analysis/Multifly analysis");
fname = sprintf('%s_tsne.pdf', db.tag);
%%
print(figure(1), fullfile(savedir, fname), '-dpdf', '-painters');
%%
for i = 2:numel(figs)
   figs(i).Position(3:4) = [560 973];    
   export_fig(figs(i), fullfile(savedir, fname), '-append', '-transparent');
end