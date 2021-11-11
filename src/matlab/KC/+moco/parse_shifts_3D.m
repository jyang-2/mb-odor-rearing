function [shifts_x, shifts_y, shifts_z] = parse_shifts_3D(shifts, Y)
%PARSE_SHIFTS_3D Summary of this function goes here
%   Detailed explanation goes here
T = length(shifts);
shifts_nr = cat(ndims(shifts(1).shifts)+1,shifts(:).shifts);
shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
shifts_x = squeeze(shifts_nr(:,2,:))';
shifts_y = squeeze(shifts_nr(:,1,:))';
shifts_z = squeeze(shifts_nr(:,3,:))';
end

