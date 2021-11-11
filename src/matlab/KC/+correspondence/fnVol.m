classdef fnVol < handle & matlab.mixin.Copyable
    %FNVOL Summary of this class goes here
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
        
        Dxw      % displaced world coordinates
        Dyw
        Dzw
        Dw
        
        V        % displacement volume
        
        A        % spatial footprint
        cm
        wcm
    end
    
    methods
        function obj = fnVol(siz,res)
            %FNVOL Construct an instance of this class
            FN = make_space_reference(siz, res);
            
            obj.siz = siz;
            obj.res = res;
            obj.R_sitk = FN.R_sitk;
            obj.xx = FN.xx;
            obj.yy = FN.yy;
            obj.zz = FN.zz;
            
            obj.xw = FN.xw;
            obj.yw = FN.yw;
            obj.zw = FN.zw;
        end
        
        function Dw = getDisplacedWorldCoordinates(obj,V)
            %GETDISPLACEDWORLDCOORDINATES get displaced world coordinates using displacement field V
            %   Detailed explanation goes here
            
            [obj.Dxw, obj.Dyw, obj.Dzw] = get_displaced_coordinates(V, obj.xw, obj.yw, obj.zw);
            Dw = [obj.Dxw(:) obj.Dyw(:) obj.Dzw(:)];
            
            obj.Dw = Dw;
            obj.V = V;

        end
        
        function wcm = getCentroidsInDisplacedWorldCoordinates(obj, A)
            %GETCENTROIDSINDISPLACEDWORLDCOORDINATES 
            %   get calculate CNMF centroids from spatial footprints, displaced 
            %   world coordinates
            wcm = get_CNMF_centroids_in_world(A, obj.Dxw, obj.Dyw, obj.Dzw);
            obj.wcm = wcm;
            obj.A = A;
        end
        
        function wcm = getComInDisplacedWorldCoordinates(obj, A)
            %GETCENTROIDSINDISPLACEDWORLDCOORDINATES 
            %   get calculate CNMF centroids from spatial footprints, displaced 
            %   world coordinates
            cm = com(A, obj.siz(1), obj.siz(2), obj.siz(3));
            %wcm = get_CNMF_centroids_in_world(A, obj.Dxw, obj.Dyw, obj.Dzw);
            wcm_x = interp3(obj.xx, obj.yy, obj.zz, obj.Dxw, cm(:,2), cm(:,1), cm(:,3)); 
            wcm_y = interp3(obj.xx, obj.yy, obj.zz, obj.Dyw, cm(:,2), cm(:,1), cm(:,3)); 
            wcm_z = interp3(obj.xx, obj.yy, obj.zz, obj.Dzw, cm(:,2), cm(:,1), cm(:,3)); 
            wcm = [wcm_x wcm_y wcm_z];
            obj.wcm = wcm;
            obj.A = A;
        end

    end
end

