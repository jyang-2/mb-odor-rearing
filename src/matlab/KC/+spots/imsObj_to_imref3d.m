function [R, imageSize, spacing] = imsObj_to_imref3d(imsObj,varargin)
%IMSOBJ_TO_IMREF3D returns an imref3D object generated from an imsObj DataSet
%                  [R, imageSize, spacing] = imsObj_to_imref3d(imsObj,varargin)
% INPUTS
%     imsObj : imaris object, like the one returned by ImarisReader('*.ims');
%         idx: [Optional, scalar], index of the DataSet 
%                when not entered, function assumes there is only one dataset
% OUTPUTS
%          R : imref3D object
%  imageSize : [xSize, ySize, zSize]
%    spacing : [xRes, yRes, zRes]


if nargin==1
   idx = 1;
elseif nargin==2
   idx = varargin{2};
end

DS = imsObj.DataSet(idx); % imsObj DataSet

%% get image size, world limits
xSize = DS.SizeX;
ySize = DS.SizeY;
zSize = DS.SizeZ;
imageSize = [xSize, ySize, zSize];

xWorldLim = [DS.ExtendMinX DS.ExtendMaxX];
yWorldLim = [DS.ExtendMinY DS.ExtendMaxY];
zWorldLim = [DS.ExtendMinZ DS.ExtendMaxZ];

%% make an imref3D object for zstack world coordinates

R = imref3d(imageSize, xWorldLim, yWorldLim, zWorldLim);

%% calculate voxel spacing/resolution

spacing = [R.PixelExtentInWorldX, R.PixelExtentInWorldY, R.PixelExtentInWorldZ];


end

