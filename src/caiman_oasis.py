import logging
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.ndimage.filters import gaussian_filter
from tifffile.tifffile import imsave

import caiman as cm
from caiman.utils.visualization import nb_view_patches3d
import caiman.source_extraction.cnmf as cnmf
from caiman.paths import caiman_datadir

from pathlib import *
# %%
import caiman as cm
from caiman.source_extraction import cnmf
from caiman.utils.utils import download_demo
from caiman.motion_correction import MotionCorrect, tile_and_correct, motion_correction_piecewise
import matplotlib.pyplot as plt

tiff_folder = Path("/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/"
                   "data/processed_data/2021-08-11/1/movie_001_ds5")
file_names = sorted(tiff_folder.glob("stk*.tif"))
for f in file_names:
    print(f)
# %%
single_movie = cm.load(str(file_names[0]))
print(single_movie.shape)

single_movie = single_movie.T
single_movie = np.moveaxis(single_movie, -1, 0)

fname = tiff_folder.joinpath('stk_001_reordered.tif')
single_movie.save(fname)


# %%
if 'dview' in locals():
    cm.stop_server(dview=dview)
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=False)
# %%
fname_new = cm.save_memmap([str(f) for f in file_names], base_name='memmap_',
                           order='C', border_to_0=0, dview=dview)
Yr, dims, T = cm.load_memmap(fname_new)
images = Yr.T.reshape((T,) + list(dims), order='F')
images = images.swapaxes(1, 3)

Cn = cm.local_correlations(images, swap_dim=False)
plt.imshow(Cn.max(0), cmap='gray')

# %%

opts_dict = {'fnames': fname,
             'strides': (48, 24, 6),  # start a new patch for pw-rigid motion correction every x pixels
             'overlaps': (24, 24, 2),  # overlap between pathes (size of patch strides+overlaps)
             'max_shifts': (4, 4, 2),  # maximum allowed rigid shifts (in pixels)
             'max_deviation_rigid': 5,  # maximum shifts deviation allowed for patch with respect to rigid shifts
             'pw_rigid': False,  # flag for performing non-rigid motion correction
             'is3D': True}

opts = cnmf.params.CNMFParams(params_dict=opts_dict)
mc = cm.motion_correction.MotionCorrect(str(fname), dview=dview, **opts.get_group('motion'))

mc.motion_correct(save_movie=True)
m_rig = cm.load(mc.fname_tot_rig, is3D=True)
# %%
mc = MotionCorrect(file_names[0], dview=dview, max_shifts=max_shifts,
                   strides=strides, overlaps=overlaps,
                   max_deviation_rigid=max_deviation_rigid,
                   shifts_opencv=shifts_opencv, nonneg_movie=True,
                   border_nan=border_nan, is3D=True)

mc.motion_correct(save_movie=True)
#%%
fname_new = "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/data/processed_data/2021-08-11/1/movie_001_ds5/memmap__d1_256_d2_256_d3_16_order_C_frames_179_.mmap"
Yr, dims, T = cm.load_memmap(str(fname_new))
images = np.reshape(Yr.T, [T] + list(dims), order='F')



# set parameters
rf = 25  # half-size of the patches in pixels. rf=25, patches are 50x50
stride = 10  # amount of overlap between the patches in pixels
K = 10  # number of neurons expected per patch
gSig = [4, 4, 2]  # expected half size of neurons
merge_thresh = 0.8  # merging threshold, max correlation allowed
p = 2  # order of the autoregressive system


cnm = cnmf.CNMF(n_processes, k=K, gSig=gSig, merge_thresh=merge_thresh, p=p, dview=dview,
                rf=rf, stride=stride, only_init_patch=False)
cnm = cnm.fit(images)