classdef tracePlot < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetObservable)
        ncids = 50;     % # of trace per figure window
        cids = 1:50;    % cell IDs to plot
        xl              % time limits (xlim)
        yl = [0 2]      % dF/F plot limits (ylim)
        axw             % width of axes (in inches)
        axh             % height of axes (in inches)
        spacing
        harray
        
    end
    
    properties        
        fig
        K = []
        T = []
        
        DFF
        hl
        hscatter
        clabels
        cidax
        
        response_ind = [];
        response_peak = [];
        stim_on
        stimlabels
        hstim
        hstimlabels
               
        axsb
        axt
        
        datestr
        fly
        nam
        htitle
        
        axprops
    end
    
    events 
        
    end
    methods
        function obj = tracePlot(DFF, varargin)
            
            addlistener(obj,'cids','PostSet',@obj.updateTraces);
            addlistener(obj, {'yl', 'xl'}, 'PostSet', @obj.updateAxLim);
            addlistener(obj, {'axw', 'axh'}, 'PostSet', @obj.updateAxSiz);
            addlistener(obj, 'spacing', 'PostSet', @obj.updateSpacing);
            
            obj.DFF = DFF;
            [K,T] = size(DFF);
            obj.K = K;
            obj.T = T; 
            obj.xl = [0 T];
            obj.axprops = {'YLim',obj.yl, 'Color', 'none', 'Clipping', 'off', 'Box', 'off',...
              'XLim', obj.xl, 'YTick', obj.yl, 'Units', 'inches'};
           
          % make figure
            obj.fig = figure('Renderer', 'painters', 'Units', 'inches');
            obj.fig.Position(3:4) = [5 9];
            
            % make trace axes
            obj.harray = gobjects(obj.ncids,1);
            obj.hl = gobjects(obj.ncids,1);
            obj.hscatter = gobjects(obj.ncids,1);
            for i = 1:obj.ncids
               ax = subplot(obj.ncids,1,i);
               obj.harray(i) = ax;
               obj.hl(i) = plot(ax, obj.DFF(obj.cids(i),:), 'k-');
               hold on
               obj.hscatter(i) = scatter(ax, [], [], 'r.');
               hold off
            end 
            
            % make harray axes invisible
            set([obj.harray(:).XAxis], 'Visible', 'off');
            set([obj.harray(:).YAxis], 'Visible', 'off');
            
            % automatic size of axes, spacing
            set(obj.harray, obj.axprops{:});
            obj.axw = obj.harray(1).Position(3);
            obj.axh = obj.harray(1).Position(4);
            obj.spacing = obj.harray(1).Position(2)-obj.harray(2).Position(2) - obj.harray(1).Position(4);

            % make scalebar axes
            obj.axsb = copyobj(obj.harray(1), obj.harray(1).Parent);
            cla(obj.axsb);
            obj.axsb.Position(2) = obj.harray(1).Position(2) + obj.spacing + obj.axh;
            obj.axsb.XAxis.Visible = 'off';
            obj.axsb.YAxis.Visible = 'on';
            obj.axsb.YLabel.String = 'dF/F';
            set(obj.axsb, obj.axprops{:});
            hold(obj.axsb);
            
            % make time scale axes
            obj.axt = copyobj(obj.harray(end), obj.harray(end).Parent);
            cla(obj.axt);
            obj.axt.XAxis.Visible = 'on';
            obj.axt.YAxis.Visible = 'off';
            obj.axt.Position(2) = obj.harray(end).Position(2)-obj.spacing - obj.axh;
            obj.axt.XLabel.String = 'time (frames)';
            set(obj.axt, obj.axprops{:});
            
            % link axt (and axsb) position to the bottom (and top) axis
            addlistener(obj.harray(end), 'Position', 'PostSet', @obj.moveAxt);
            addlistener(obj.harray(1), 'Position', 'PostSet', @obj.moveAxsb);
            
            obj.cidax =  gobjects(obj.ncids,1);
            obj.clabels =  gobjects(obj.ncids,1);
            for i = 1:obj.ncids
                pos = obj.harray(i).Position;
                pos(1) = pos(1) - 0.05;
                pos(3) = 0;
                obj.cidax(i) =  axes(obj.fig, 'Units', 'inches', 'Position', pos,...
                    'Color', 'none', 'XTick', [], 'YTick', [], 'Visible', 'off');
                obj.cidax(i).XAxis.Visible = 'off';
                obj.cidax(i).YAxis.Visible = 'off';
                linkaxes([obj.cidax(i) obj.harray(i)], 'y');
                obj.clabels(i) = text(obj.cidax(i), 0, 0, num2str(i), 'FontSize', 7, ...
                    'HorizontalAlignment', 'right',  'Visible', 'off');
            end
            set(obj.cidax, 'Visible', 'on');
            set(obj.clabels, 'Visible', 'on');
        end
        
        function updateTraces(obj, src, event)
           for i = 1:numel(obj.cids)
               obj.hl(i).YData = obj.DFF(obj.cids(i),:);
               obj.clabels(i).String = num2str(obj.cids(i));
           end
           obj.updateScatter();
           if numel(obj.cids) < obj.ncids
              for i = numel(obj.cids)+1:obj.ncids
               obj.hl(i).YData = [];
               obj.clabels(i).String = '';
              end
           end
        end
        
        function setTitle(obj,tt)
            obj.htitle = sgtitle(tt,'FontSize', 10);
            drawnow;
        end
        
        function showStimOn(obj,src,event)
            if ~isempty(obj.hstim)
                delete(obj.hstim);
            end
            
            if ~isempty(obj.stim_on)
                x = repelem(obj.stim_on(:), 1,2);
                y = [0 -(obj.spacing+obj.axh)*obj.yl(2)*(obj.ncids+2)/(obj.axh)];
                y = repmat(y, size(x,1),1);
                disp([x y])
                obj.hstim = plot(obj.axsb, x',y', 'r:');
            end
            drawnow;
        end
        
        function showStimLabels(obj, src, event)
            if ~isempty(obj.hstimlabels)
                delete(obj.hstimlabels);
            end
            
            if ~isempty(obj.stimlabels)
                axes(obj.axsb);
                obj.hstimlabels = text(obj.stim_on, ones(size(obj.stim_on)), obj.stimlabels, ...
                    'FontSize', 7, 'Rotation', 45);
            end
            drawnow;
        end
        
        function updateSpacing(obj, src, event)
           for i = obj.ncids-1:-1:1
               obj.harray(i).Position(2) = obj.harray(i+1).Position(2) + obj.spacing + obj.harray(i).Position(4);
               obj.cidax(i).Position(2) = obj.harray(i+1).Position(2) + obj.spacing + obj.harray(i).Position(4);
           end
           obj.showStimOn;
           obj.showStimLabels;
           drawnow
         end
         
        
        function updateAxLim(obj,src,event)
            for i = 1:obj.ncids
                obj.harray(i).XLim =  obj.xl;
                obj.harray(i).YLim =  obj.yl;
               
            end
            obj.axsb.XLim = obj.xl;
            obj.axt.XLim = obj.xl;
            obj.axsb.YLim = obj.yl;
            obj.axt.YLim = obj.yl;
            obj.axsb.YTick = obj.yl;
            obj.showStimOn;
            obj.showStimLabels;
            drawnow
        end        
       
        function moveAxt(obj,src,event)
            drawnow
            disp('dostuff');
            obj.axt.Position(2) = obj.harray(end).Position(2)-2*obj.axh;
            drawnow
        end
        
        function moveAxsb(obj,src,event)
            drawnow
            disp('dostuff');
            obj.axsb.Position(2) = obj.harray(1).Position(2) + 2*obj.axh;
            drawnow
        end
         
        function updateAxSiz(obj,src,event)
            snam = src.Name;
            switch snam
                case 'axh'
                    for i = 1:obj.ncids
                        obj.harray(i).Position(4) = obj.axh;
                        obj.axsb.Position(4) = obj.axh;
                        obj.axt.Position(4) = obj.axh;
                    end

                case 'axw'
                    for i = 1:obj.ncids
                        obj.harray(i).Position(3) = obj.axw;
                        obj.axsb.Position(3) = obj.axw;
                        obj.axt.Position(3) = obj.axw;
                    end
            end
            drawnow
        end
        
        function setResponseInfo(obj, response_ind, response_peak)
            obj.response_ind = response_ind;
            obj.response_peak = response_peak;
        end
        
        function updateScatter(obj, src, event)
           for i = 1:numel(obj.cids)
               cid = obj.cids(i);
               set(obj.hscatter(i), 'XData', obj.response_ind(cid,:), 'YData', obj.response_peak(cid,:)+0.4);
           end
           for i = (numel(obj.cids)+1):obj.ncids
               set(obj.hscatter(i), 'XData', [], 'YData', []);
           end
        end
        
        function bottomSpacing(obj, x)
            %allpos=cell2mat(get(findobj('-property', 'TightInset'), 'TightInset'));
            %dp = min(allpos(:,2))-x;
            dp = obj.axt.Position(2)-x;
            
            for i = 1:obj.ncids
                obj.harray(i).Position(2) = obj.harray(i).Position(2) - dp;
                obj.cidax(i).Position(2) = obj.cidax(i).Position(2) - dp;
            end
            %obj.axsb.Position(2) = obj.axsb.Position(2) - dp;
            %obj.axt.Position(2) = obj.axt.Position(2) - dp;
        end
        function leftSpacing(obj, x)
            dp = obj.axsb.Position(1)-x;
            for i = 1:obj.ncids
                obj.harray(i).Position(1) = obj.harray(i).Position(1) - dp;
                obj.cidax(i).Position(1) = obj.cidax(i).Position(1) - dp;

            end
            obj.axsb.Position(1) = obj.axsb.Position(1) - dp;
            obj.axt.Position(1) = obj.axt.Position(1) - dp;
        end
        
        
    end
end

