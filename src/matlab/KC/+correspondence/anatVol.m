classdef anatVol < handle & matlab.mixin.Copyable
    %ANATVOL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        siz      % size in px
        res      % px resolution
        R_sitk   % imref3D for simpleITK 
        
        xx       % meshgrid, intrinsix x-coords
        yy       % meshgrid, intrinsix y-coords
        zz       % meshgrid, intrinsix z-coords
        
        xw       % meshgrid, world x-coords
        yw       % meshgrid, world y-coords
        zw       % meshgrid, world z-coords
        
        posXYZ
        cm_px
        wcm
        labels
        seeds
        
        imsObj 
        imsInfo
        R_ims
    end
    
    methods
       function obj = anatVol(siz,res)
            %FNVOL Construct an instance of this class
            S = make_space_reference(siz, res);
            
            obj.siz = siz;
            obj.res = res;
            obj.R_sitk = S.R_sitk;
            
            obj.xx = S.xx;
            obj.yy = S.yy;
            obj.zz = S.zz;
            
            obj.xw = S.xw;
            obj.yw = S.yw;
            obj.zw = S.zw;
        end
        
        function [wcm, cm] = loadImarisSpots(obj,imsObj, spot_idx)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.posXYZ = imsObj.Spots(1, spot_idx).GetPositions();
            obj.imsInfo = spots.read_ims_dataset(imsObj);
            [obj.R_ims,~,~] = spots.imsObj_to_imref3d(imsObj);
            [obj.seeds, obj.labels, cm] = spots.make_spot_seeds(obj.posXYZ, obj.R_ims);   

            [xI, yI, zI] = worldToIntrinsic(obj.R_ims, obj.posXYZ(:,1), obj.posXYZ(:,2), obj.posXYZ(:,3));
            [xW, yW, zW] = intrinsicToWorld(obj.R_sitk, xI,yI,zI);
            
            wcm = [xW, yW, zW];
    
            obj.cm_px = cm;
            obj.wcm = wcm;
        end
    end
end

