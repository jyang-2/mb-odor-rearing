classdef LabelmatDisplay < handle & matlab.mixin.SetGet
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        L
        I
        labelstats
        labels
        cm
        label_planes
        ax
        h_im
        h_fig
        
        siz
        
        active_label
        active_plane

        clim
    end
    
    methods
        function obj = LabelmatDisplay(I, L, ax)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.I = I;
            obj.L = L;
            obj.siz = size(I);

            if nargin>2
                if isgraphics(ax)
                   obj.ax = ax; 
                end
            end
            obj.clim = [1400 2900];
            obj.labels = unique(L(:));
            obj.labels(obj.labels==0)=[];
            obj.labelstats = regionprops3(L);
            obj.cm = obj.labelstats.Centroid;
            obj.cm = obj.cm(:,[2 1 3]);
            obj.label_planes = round(obj.cm(:,3));
            
        end
        
        function initializeDisplay(obj)
            %SETAXISHANDLE assigns axis handle to draw
            %   Detailed explanation goes here
                       % draw background image

           if isempty(obj.ax)
              figure('Name', 'Labelmatrix display');
              set(obj, 'ax', gca);
              obj.h_fig = ancestor(obj.ax, 'figure');
           end
           
           if ~isgraphics(obj.ax)
             figure('Name', 'Labelmatrix display');
              set(obj, 'ax', gca);
              obj.h_fig = ancestor(obj.ax, 'figure');

           end
           
           obj.active_plane = 1;
           obj.h_im = imagesc(obj.ax, obj.I(:,:,obj.active_plane), obj.clim);          
           colormap gray; axis equal; axis tight; 
           %axis off;

           set(gcf, 'WindowButtonDownFcn', @obj.mouseClick);  % mouse clicks
              set(gcf, 'WindowScrollWheelFcn', @obj.scrollWheel);  % mouse clicks

        end

        function scrollWheel(obj, src, eventData)
            modifiers = get(obj.h_fig,'CurrentModifier');
            %disp(get(obj.h_fig,'CurrentModifier'))
                        
            wasShiftPressed = ismember('shift',   modifiers);  % true/false
            %wasCtrlPressed  = ismember('control', modifiers);  % true/false
            %wasAltPressed   = ismember('alt',     modifiers);  % true/false
            if wasShiftPressed
                dz = eventData.VerticalScrollCount;
                new_active_plane = obj.active_plane + dz;
                if dz<0
                    obj.active_plane = max(1,new_active_plane);
                elseif dz>0
                    obj.active_plane = min(obj.siz(end), new_active_plane);
                end
                if ~isempty(obj.active_label)
                    obj.drawLabel(obj.active_label, obj.active_plane);
                else
                    set(obj.h_im, 'CData', obj.I(:,:,obj.active_plane));
                end
                    
            end
        end
        

        function mouseClick(obj, src, eventData)
            modifiers = get(obj.h_fig,'CurrentModifier');
            %disp(get(obj.h_fig,'CurrentModifier'))
                        
            wasShiftPressed = ismember('shift',   modifiers);  % true/false
            wasCtrlPressed  = ismember('control', modifiers);  % true/false
            wasAltPressed   = ismember('alt',     modifiers);  % true/false
            
            if wasShiftPressed
                [x, y] = ginput(1);
                hold on
                %plot(obj.ax, x(1), y(1), 'ro')
                %drawnow
                d = pdist2([y,x,obj.active_plane], obj.cm);
                new_label = find(d==min(d));
                disp(new_label);
                obj.drawLabel(new_label);
            end
        end

        function drawLabel(obj,ll,z)
            if ~exist('z', 'var')
                z = obj.label_planes(ll);
                obj.active_plane = z;
            end
            
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            %z = obj.label_planes(ll);
            im = mat2gray(obj.I(:,:,z), obj.clim);
            imoverlay = labeloverlay(im, obj.L(:,:,z), 'IncludedLabels', ll, 'Transparency', 0.4, 'Colormap', 'parula');
            %imshow(obj.ax,imoverlay);
            set(obj.h_im, 'CData', imoverlay);
           colormap gray; axis equal; axis tight; 
           %axis off;
           obj.active_label = ll;

            
            str = sprintf('Label # %d (z=%d)', ll, z);
            
            title(obj.h_im.Parent, str, 'FontWeight', 'normal', 'FontSize', 10);
        end
    end
end

