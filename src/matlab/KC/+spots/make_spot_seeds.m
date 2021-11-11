function [seeds, labels, cm_px] = make_spot_seeds(posXYZ, R)
%MAKE_SPOT_SEEDS    takes imaris spot coordinates and imref3D object, and
%                   returns a seed labelmatrix, a list of labels, and
%                   centroids in pixel coordinates.
%
%                   [seeds, labels, cm_px] = make_spot_seeds(posXYZ, R)
% 
% INPUTS
%     posXYZ : imaris spot coordinates (like Spots(1).GetPositions()
%          R : imref3D object (create it using the imsObj_to_imref3d fcn)
%
% OUTPUTS
%      seeds : [m x n x k] labelmatrix, w/ seeds of spot centers
%     labels : list of spot labels in seeds
%      cm_px : [n x 3] array of spot center coordinates



[i,j,k] = worldToSubscript(R, posXYZ(:,1), posXYZ(:,2), posXYZ(:,3));
cm_px  = [i,j,k];
numSpots = size(posXYZ,1);

seeds = zeros(R.ImageSize);
for ind =1:numSpots
    if ~isnan(i(ind)) && ~isnan(j(ind)) && ~isnan(k(ind))
        seeds(i(ind), j(ind), k(ind))=ind; 
    end
end

labels = unique(seeds);
labels(labels==0) = [];

seeds = uint16(seeds);
end