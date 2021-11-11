% generate a label seed image from an imaris file

%%

xPath = 'demo_zstack.ims';
imsObj = ImarisReader(xPath);
ims_info = spots.read_ims_dataset(imsObj);

%% get imref3D, image size, and image resolution
[R,siz,res] = spotsimsObj_to_imref3d(imsObj);
%%
posXYZ = imsObj.Spots(1).GetPositions;
[seeds, labels, cm_px] = spots.make_spot_seeds(posXYZ, R);


