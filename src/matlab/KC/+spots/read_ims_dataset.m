function ims_info = read_ims_dataset(imsObj)
%READ_IMS_DATASET returns info about imsObj
%
% INPUTS
%     imsObj : imaris object, like the one returned by ImarisReader('*.ims');
%
% OUTPUTS
%   ims_info : struct, contains information about contents of imsObj
%%
ims_info = struct;
p = string(properties(imsObj.DataSet))';
for str = p
   %disp(str)
   ims_info.(str) = get(imsObj.DataSet, str);
   %disp(ims_info.(str));
   
end

end

