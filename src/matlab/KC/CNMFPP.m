classdef CNMFPP
    %CNMFPP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = CNMFPP(inputArg1,inputArg2)
            %CNMFPP Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    methods(Static)
        %%
        function cids = getCenterIDs(centers, pts)
            %getCenterIDs....returns the indices of the corresponding
            %components given a list of points (point = cm (center))
            
           [~, ia, ~] = intersect(centers, pts);
           cids = ia;         
        end
        
        function cids = getComponentIDs(A, siz, pts)
            centers = com(A, siz(1),siz(2),siz(3));
            
           [~, ia, ~] = intersect(centers, pts);
           cids = ia;       
        end
        
        function [cids, bcids] = getMaskedCenterIDs(centers, mask, label)
           mask = (mask==label);
           centers = round(centers);
           
           idx = sub2ind(size(mask), centers(:,1), centers(:,2), centers(:,3));
           
           bcids = mask(idx);
           cids = find(bcids);
        end
        
        function normedA = getNormedSpatialComponents(A)
            %returns an L2-normed A (spatial components of CNMF)
            nA = sqrt(sum(A.^2,1));
            normedA = bsxfun(@times, A, 1./nA);
        end
        
        function spatial = getSpatialImage(A,siz)
            spatial = full(sum(A,2));
            spatial = reshape(spatial, siz);
        end
        
        function spatial = getNormedSpatialImage(A, siz)
            A = CNMFPP.getNormedSpatialComponents(A);
            spatial = CNMFPP.getSpatialImage(A, siz);
            
        end
        
        function Hmw = plotTraces(traces, cids, ti, response, title_text)
           % traces with responders marked (if response vectors were calculated)
           
            T = size(traces,2);
            figure
            clf;
            Hmw = iosr.figures.multiwaveplot(1:T, [], traces(cids,:),...
                'reverseY', true,...
                'gain', 5,...
                'mode', 'plot');
            hold on
            
            %if ~exist('options', 'var') || ~isfield('options', 'tt')
            
            if exist('ti', 'var') && ~isempty(ti)
               ax = gca;
                ax.YTick = 1:numel(cids);
                ax.YTickLabel = cids;
                ax.FontSize=8;
                ax.TickDir = 'out';
                ax.YLabel.String = 'dF/F';
                %ax.Box = 'off';

                set(gca, 'XTick', ti.stim_on,...
                    'XTickLabel', ti.pin_odors(ti.si),...
                    'XAxisLocation', 'top');    
                ax.YAxis.FontSize=8;
                ax.XAxis.FontSize=8;
                xtickangle(90)
            end

            if exist('response', 'var') && ~isempty(response)
                vline(ti.stim_on);
                %vline(ti.peakresp_start, 'k--')
                vline(ti.peakresp_end)
                [r,c] = find(response.bin_response(cids,:));
                x= nan(size(r));
                y = nan(size(c));
                for i = 1:numel(r)
                    ri = r(i);
                    ci = response.peak_ind(cids(ri),c(i));
                    x(i) = Hmw(ri).XData(ci);
                    y(i) = Hmw(ri).YData(ci);
                end
                scatter(x,y-0.15, 10,'ro', 'filled');

            end

            if exist('title_text', 'var') && ~isempty(title_text)
                title(title_text, 'Interpreter', 'none') 
            end
            
            
        end
        
        function fig = plotOverlaidTraces(traces1, traces2, cids, ti, response, title_text)
            T = size(traces1,2);

            fig = figure('Visible', 'on');
            set(gcf,'renderer', 'painters')
            colormap(hsv)

            LL = gobjects(numel(cids),1);
            
            spacing = 3;
            scaling = 1;
            for i = 1:numel(cids)
                hold on
                text(-25, -i*spacing, num2str(cids(i)), 'FontSize', 7, 'HorizontalAlignment', 'right');
                l1 = plot((1:T), traces1(cids(i),:)*scaling-i*spacing, 'Color', 'k', 'LineWidth', 1);
                l1.Color(4) = .3;
                LL(i) = l1;
                plot((1:T), traces2(cids(i),:)*scaling-i*spacing, 'Color', 'k', 'LineWidth', 1);
            end
            

            set(gca, 'TickDir', 'out', 'FontSize', 8,... 'YLabel.String', 'dF/F', ...
                'YTick', []);

           
             if exist('ti', 'var') && ~isempty(ti)
               ax = gca;
                ax.FontSize=8;
                ax.TickDir = 'out';
                ax.YLabel.String = 'dF/F';
                %ax.Box = 'off';

                set(gca, 'XTick', ti.stim_on,...
                    'XTickLabel', ti.pin_odors(ti.si),...
                    'XAxisLocation', 'top');    
                ax.YAxis.FontSize=8;
                ax.XAxis.FontSize=8;
                ax.YTick = [];
               
                xtickangle(90)
             end
            
             if exist('response', 'var') && ~isempty(response)
                vline(ti.stim_on);
                [r,c] = find(response.bin_response(cids,:));
                x= nan(size(r));
                y = nan(size(c));
                for i = 1:numel(r)
                    ri = r(i);
                    ci = response.peak_ind(cids(ri),c(i));
                    x(i) = LL(ri).XData(ci);
                    y(i) = LL(ri).YData(ci);
                end
                scatter(x,y+0.5, 10,'ro', 'filled');

             end
            
             if exist('title_text', 'var') && ~isempty(title_text)
                title(title_text, 'Interpreter', 'none') 
             end
             
             x0 = LL(1).XData(1);
             y0 = LL(1).YData(1);
             sb = line([x0 x0], [y0+2 y0+2+2*scaling], 'LineWidth', 1, 'Color', 'r');
             text(x0+10, y0+2+scaling, '200%', 'FontSize', 7);
                         
             ax = gca;
             ax.XLim = [-0.08 1.05] * T;
             ax.YLim = [min(LL(end).YData)-2, max(LL(1).YData)+5];
        end
    end
end

