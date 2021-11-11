function run_cnmf_in_folder(folder)
%RUN_CNMF_IN_FOLDER Summary of this function goes here
%   Detailed explanation goes here
cnmf_dir = fullfile(folder, 'cnmf');
fname_opts = fullfile(cnmf_dir, 'cnmfoptions.mat');
cnmfoptions = load(fullfile(cnmf_dir, 'cnmfoptions.mat'));


tiff_listing = dir(fullfile(folder, '*.tif'));
assert(length(tiff_listing)==1, sprintf('The folder (%s) contains more than 1 tif file', folder));

Y = loadtiff(fullfile(tiff_listing.folder, tiff_listing.name));

CNM = CNMF();
CNM.loadArray(Y);
CNM.optionsSet(cnmfoptions.options);                        % setup the options structure
CNM.gSig = cnmfoptions.options.gSig;

CNM.preprocess();             % preprocessing (compute some quantities)
CNM.initComponents(cnmfoptions.K);      % initialization
%CNM.plotCenters()           % plot center of ROIs detected during initialization

CNM.updateSpatial();        % update spatial components
CNM.updateTemporal(0);      % update temporal components (do not deconvolve at this point)

%===============plot contours===================
CNM.cm = [];
CNM.contours=[];
CNM.COM();
CNM.K = size(CNM.A,2);
mY = mean(CNM.Y,3);
CNM.CI=mY;

fig = figure();
CNM.plotContours(1);
fig.Position(3:4)=[800 800];
F = getframe(gca);
imwrite(F.cdata, fullfile(cnmf_dir, 'contours.png'));
close(fig)

save(fullfile(cnmf_dir, 'CNM.mat'), 'CNM');
end

