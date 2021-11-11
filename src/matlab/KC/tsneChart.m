classdef tsneChart < matlab.graphics.chartcontainer.ChartContainer & matlab.graphics.chartcontainer.mixin.Legend
%TSNECHART Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        XData (1,:) double = nan
        YData (1,:) double = nan
        group (1,:) 
        TitleText (:,:) char = '';
        GscatterArray
        Rois (1,:) = gobjects(0);
        newRoi (1,1) images.roi.Freehand
        selectedPts
    end
    
    properties(Access = private)
        
    end
    
    methods
        function pts = getPointsInRoi(obj, n)
            pts = checkPoints(obj.XData, obj.YData, obj.Rois(n));
        end
    end
        
    methods(Access=protected)
        function setup(obj)
           ax = getAxes(obj);
           ax.Color = 'none';
           obj.GscatterArray = gscatter(ax, obj.XData, obj.YData, obj.group); 
           %obj.newRoi = images.roi.Freehand;
           %obj.Rois = images.roi.Freehand(ax, 'LineWidth', 2, 'Label', num2str(1));
           set(obj.Parent, 'Renderer', 'painters', 'Name', 'tsne', 'tag', 'tsne')

        end
        
        function update(obj)
            ax = getAxes(obj);
            obj.Rois = findobj(ax, 'Type', 'images.roi.Freehand');
            obj.Rois = fliplr(obj.Rois);
            for i = 1:numel(obj.Rois)
               obj.Rois(i).Label = num2str(i);
            end
            obj.selectedPts = checkPoints(obj.XData, obj.YData, obj.Rois);
            delete(obj.GscatterArray )
            obj.GscatterArray = gscatter(ax, obj.XData, obj.YData, obj.group, [], [], 8, 'off');
            title(ax, obj.TitleText);
            lgd = getLegend(obj);
            %lgd.String = num2str(unique(obj.group)');
            lgd.String = string(unique(obj.group)');
            lgd.Location = 'northeastoutside';
            lgd.Color = 'none';
            lgd.Box = 'on';

        end
        

    end
    methods
        function drawRoi(obj)
           ax = getAxes(obj);
           obj.newRoi = drawfreehand(ax);
           hold(ax, 'on')
           
        end
        
        function deleteRois(obj, cids)
            %ax = getAxes(obj);
            delete(obj.Rois(cids));
            obj.Rois(cids) = [];
        end
        
        function figs = plotRoiTraces(obj, DFF, sDFF, ti, response)
            nroi = size(obj.selectedPts,2);
            figs = gobjects(0);
            traces_per_fig = 100;
            dar = [1.556010000000000e+03,81.598072978261230,1];
            
            for i = 1:nroi
                cids = find(obj.selectedPts(:,i));
                while(numel(cids)>100)
                    ft = CNMFPP.plotOverlaidTraces(DFF, sDFF, cids(1:traces_per_fig), ti, response);
                    
                    title(i);
                    cids(1:traces_per_fig)=[];
                    figs(numel(figs)+1) = ft;
                end
                ft = CNMFPP.plotOverlaidTraces(DFF, sDFF, cids, ti, response);
                title(i);
                figs(numel(figs)+1) = ft;
            end
            figs(i).Position(3:4) = [560 973];
            
            figs = setFigSize(figs, 560, 973);
            figs = setDAR(figs, dar);
            set(figs, 'Name', 'traces');
            
        end

    end
    
end

  
function pts_selected = checkPoints(x,y,rois)
        pts_selected = false(numel(x), numel(rois));
        for i = 1:numel(rois)
            pts_selected(:,i) = inROI(rois(i), x, y);
        end     
end

function hfigs = setFigSize(hfigs, w,h)
    for i = 1:numel(hfigs)
       hfigs(i).Position(3:4) = [w h];
    end
end

function hfigs = setDAR(hfigs, rat)
    hax = findobj(hfigs, 'Type', 'axes');
    set(hax, 'DataAspectRatio', rat);
end

function printPdfs(obj, hfigs, savedir)
    fname = sprintf('%s_tsneChart.pdf')
    print(figure(1), fullfile(savedir, fname), '-dpdf', '-painters');

end
