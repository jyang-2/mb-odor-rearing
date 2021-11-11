A = CNMc.A;
siz = CNMc.dims; 
[d1,d2,d3] = matsplit(siz);
CNMc.COM();
cm = CNMc.cm;
Cn = reshape(CNMc.sn, CNMc.dims);


%%
K = size(A,2);
best_planes = get_component_plane(A, siz);

m = floor(sqrt(d3));
n = ceil(sqrt(d3));
XData = cm(:,2);
YData = cm(:,1);
ZData = cm(:,3);


%tl = tiledlayout(m,n);
tl = tiledlayout('flow');
tl.TileSpacing = 'none';
tl.Padding = 'compact';

% draw background images by plane
harray = gobjects(1,z);
for z = 1:d3
   nexttile;
   harray(z) = gca;
   imagesc(Cn(:,:,z));
   colormap gray
   axis image
   title(num2str(z));
   hold on
   custom_roi_toolbar(harray(z));
end

roilist = gobjects(1,K);
for i = 1:K
    z = best_planes(i);
    roilist(i) = plot(harray(z), XData(i), YData(i), '.');
end
%%
cids_selected = false(K,1);
for i = 1:d3
    fhs = findobj(harray(i), 'Type', 'images.roi.rectangle');
    
    idx = find(best_planes==i);
    pts_selected = false(size(idx));
    for ir = 1:numel(fhs)
        fh = fhs(ir);
        pts_selected = pts_selected | inROI(fh, XData(idx), YData(idx));
    end
    cids_selected(idx) = pts_selected;
end
%%
keep = true(size(CNMc.A,2),1);
keep(cids_selected) = false;

CNMc.keepComponents(keep);