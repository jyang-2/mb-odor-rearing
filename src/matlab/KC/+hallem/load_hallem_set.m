function hallem_set = load_hallem_set()
%   LOAD_HALLEM_SET     Returns struct "hallem"
%  Outputs:
%  hallem
%   receptor
% 		all: [1x24 string], names of all glomeruli
% 		phr: [1x4 string], names of pheremone glomeruli
%     nonphr: [1x20 string], names of non-pheremone glomeruli% 
%
%   chem
% 	    cat_list: [10x1 string], chemical categories
% 	      	cmap: [10x3 double], colormap of chem. categories, out of 255
% 	    cmap_rgb: [10x3 double], colormap of chem. categories, out of 1
%
% 	odors
% 	   odor_list: [111x1 string], names of hallem odors
% 	      abbrev: [111x1 string], abbreviations for hallem odors

%p = mfilename('fullpath');
%[xpath, ~] = fileparts(mfilename('fullpath'));
%hallem_set = load(fullfile(xpath, 'hallem.mat'));

hallem_set = load('hallem-8ea3408c-de8c-42a6-a609-8bf4c3770264.mat');
cmaps = load('favorite_colormaps.mat');
hallem_set.cmaps = cmaps;
end

