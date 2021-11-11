cmap = [166,206,227;
31,120,180;
178,223,138;
51,160,44;
251,154,153;
227,26,28;
253,191,111;
255,127,0;
202,178,214];
cmap = cmap/255;
%%

clr_stack = ones(256,256,length(cmap)) .* reshape([1:length(cmap)], 1, 1, length(cmap));

%%
cgrp = 0;
normedA = CNMFPP.getNormedSpatialComponents(CNMc.A(:,cnodes==cgrp));
%%
spatial = CNMFPP.getSpatialImage(normedA,[256 256 9]);

%%
z = 1;
h2 = axes;
%%
h_ov = imagesc(clr_stack(:,:,z));
hold off
set(h_ov, 'AlphaData', spatial(:,:,z));

           
  
%%           
set(obj.h_ov, 'CData', rgb_i, 'AlphaData', alpha_map);
