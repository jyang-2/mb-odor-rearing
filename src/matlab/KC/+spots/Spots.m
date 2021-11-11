classdef Spots < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        num_spots
        labels   % labels to be used
        spot_stats % spot statistics , includes intensity info
        cm_pos
        cm_sub
        cm_ind
        
        ims_info
        exp_info
        
        spacing
        siz
        dim
        
        extentX    %extent in world coordinates
        extentY
        extentZ
        
        R % imref3d object
    end
    
    methods
        function obj = Spots()
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
        end
        
%         function readImsInfo(obj, imsObj)
%             obj.ims_info = read_ims_dataset(imsObj);
%         end
        
        function readExpInfo(obj,file)
           obj.exp_info = readThorImageXML(file); 
        end
        
        function readSpots(obj, imsSpots)
            obj.num_spots = imsSpots.NumberOfSpots;
            obj.labels = imsSpots.GetIDs;
            obj.spot_stats = imsSpots.GetStatistics;
            obj.cm_pos = imsSpots.GetPositions;
            obj.R = imref3d(obj.dim, obj.extentX, obj.extentY, obj.extentZ);
            [i,j,k] = worldToSubscript(obj.R, obj.cm_pos(:,1), obj.cm_pos(:,2), obj.cm_pos(:,3));
            obj.cm_sub = [i,j,k];
            obj.cm_ind = sub2ind(obj.dim, i,j,k);
        end
        
        function setSpacing(obj,S)
           if isnumeric(S)
               obj.spacing = S;
           elseif strcmp(class(S), 'ImarisReader') %#ok<STISA>
               obj.ims_info = read_ims_dataset(S);
               obj.extentX = [obj.ims_info.ExtendMinX, obj.ims_info.ExtendMaxX];
               obj.extentY = [obj.ims_info.ExtendMinY, obj.ims_info.ExtendMaxY];
               obj.extentZ = [obj.ims_info.ExtendMinZ, obj.ims_info.ExtendMaxZ];
               obj.siz = [diff(obj.extentX), diff(obj.extentY), diff(obj.extentZ)];
               obj.dim = [obj.ims_info.SizeX, obj.ims_info.SizeY, obj.ims_info.SizeZ];
               obj.spacing = obj.siz./obj.dim;
           elseif (isstring(S)||ischar(S)) && strcmp(S(end-3:end),'.xml')
               obj.exp_info=readThorImageXML(exp_info_filepath);
           end        
        end
        
        function seeds = getSeedArray(obj, labelsToUse)
            seeds = uint16(zeros(obj.dim));
            if ~exist('labelsToUse', 'var')
               seeds(obj.cm_ind) =  obj.labels;
            else
                seeds(obj.cm_ind(obj.labels==labelsToUse))=labelsToUse;
            end
        end
    end
end

