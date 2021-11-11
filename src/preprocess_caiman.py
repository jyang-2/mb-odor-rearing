#import bokeh.plotting as bpl
import cv2
import glob
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from ryim import thorimagels
from tifffile.tifffile import imsave
import caiman as cm
import caiman.source_extraction.cnmf as cnmf

#%%
flydate = '2021-08-30'
#flynum = 2
#movnam = 'movie_001'

DATA_2P_DIR = Path.cwd().joinpath('data', '2p')
xml_file = Path.cwd().joinpath('data', '2p', flydate, str(flynum), movnam, 'Experiment.xml')
raw_file = Path.cwd().joinpath('data', '2p', flydate, str(flynum), movnam, 'Image_001_001.raw')

save_dir = raw_file.parent.relative_to(DATA_2P_DIR)
save_dir = Path.cwd().joinpath('data', 'processed_data', save_dir)
print(save_dir)

meta = thorimagels.read_metadata(xml_file)

meta_dict = dict(siz=thorimagels.movie_shape(meta), )

siz = thorimagels.movie_shape(meta)

#%%
Y = thorimagels.read_raw(raw_file).reshape(siz)
Y = Y[:, :16, :, :]

Y = Y.T
Y = np.transpose(Y, (-1, 0, 1, 2))
# %%
T, d1, d2, d3 = Y.shape
dims = (d1, d2, d3)
imsave(save_dir.joinpath(fname), )

# %%
for i, y in enumerate(np.split(Y, 3, axis=0)):
    fname = f"stk_{str(i).zfill(4)}.tif"

    print(save_dir.joinpath(fname).relative_to(Path.cwd()))
    print(y.shape)
    imsave(save_dir.joinpath(fname).relative_to(Path.cwd()), y)
    print('saved!\n')
# %%

files = sorted(save_dir.glob("stk_*.tif"))
fnames = [str(item) for item in files]
opts_dict = {'fnames': fnames,
             'strides': (48, 48, 4),  # start a new patch for pw-rigid motion correction every x pixels
             'overlaps': (16, 16, 2),  # overlap between pathes (size of patch strides+overlaps)
             'max_shifts': (8, 8, 4),  # maximum allowed rigid shifts (in pixels)
             'max_deviation_rigid': 4,  # maximum shifts deviation allowed for patch with respect to rigid shifts
             'pw_rigid': False,  # flag for performing non-rigid motion correction
             'is3D': True,
             'num_frames_split': 323,
             'use_cuda': True}

opts = cnmf.params.CNMFParams(params_dict=opts_dict)
print(opts)
#%% set up cluster
if 'dview' in locals():
    cm.stop_server(dview=dview)
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=False)
#%% restart cluster to clean up memory
cm.stop_server(dview=dview)
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=False)
#%% run motion correction, mmap movie
mc = cm.motion_correction.MotionCorrect(fnames, dview=dview, **opts.get_group('motion'))
mc.motion_correct(save_movie=True)
# n_planes = thorimagels.z_steps(meta)
# for i in range(n_planes):
#     y = Y[:, i, :, :]
#     fname = f"plane_{str(i).zfill(4)}.tif"
#     imsave(raw_file.parent.joinpath(fname), y)

#%% plot shifts, visualize results
shift_x, shift_y, shift_z = zip(*mc.shifts_rig)
plt.plot(shift_x, label='x')
plt.plot(shift_y, label='y')
plt.plot(shift_z, label='z')
plt.legend()
plt.show()

# load mmap rigid
Y_rig = cm.load(mc.fname_tot_rig, is3D=True)

# plot template image
plt.imshow(mc.total_template_rig, cmap = 'gray')

# correlation image
Cn = cm.local_correlations(Y_rig, swap_dim=False)
#%% Piecewise rigid
mc.pw_rigid = True  # turn the flag to True for pw-rigid motion correction
mc.template = mc.mmap_file  # use the template obtained before to save in computation (optional)

mc.motion_correct(save_movie=True, template=mc.total_template_rig)
m_els = cm.load(mc.fname_tot_els)
#%%
#plt.subplot(1,2,1); plt.imshow(m_orig.local_correlations(eight_neighbours=True, swap_dim=False))
Cn_rig = cm.local_correlations(Y_rig, eight_neighbours=True, swap_dim=False)
Cn_els = cm.local_correlations(m_els, eight_neighbours=True, swap_dim=False)
plt.figure(figsize = (20,10))
plt.subplot(1,2,1); plt.imshow(Cn_rig[..., 7])
plt.subplot(1,2,2); plt.imshow(Cn_els[..., 7])
plt.show()
#%%
fname_new = cm.save_memmap(mc.fname_tot_rig, base_name='memmap_', order='C',
                           border_to_0=0, dview=dview)  # exclude borders

# now load the file
Yr, dims, T = cm.load_memmap(fname_new)
images = np.reshape(Yr.T, [T] + list(dims), order='F')
#%%
rf = 25  # half-size of the patches in pixels. rf=25, patches are 50x50
stride = 10  # amount of overlap between the patches in pixels
K = 150  # number of neurons expected per patch
gSig = [2, 2, 2]  # expected half size of neurons
merge_thresh = 0.8  # merging threshold, max correlation allowed
p = 2  # order of the autoregressive system
ssub = 1
tsub = 3
min_corr = 0.7

cnm = cnmf.CNMF(n_processes, k=K, gSig=gSig, merge_thresh=merge_thresh, p=p, dview=dview,
                rf=rf, stride=stride, only_init_patch=True, ssub=ssub, tsub=tsub, min_corr=min_corr)
cnm = cnm.fit(images)

print(('Number of components:' + str(cnm.estimates.A.shape[-1])))
#%%

fr = 3 # approx final rate  (after eventual downsampling )
decay_time = 10.  # length of typical transient in seconds
use_cnn = False  # CNN classifier is designed for 2d (real) data
min_SNR = 1.5      # accept components with that peak-SNR or higher
rval_thr = 0.6   # accept components with space correlation threshold or higher
cnm.params.change_params(params_dict={'fr': fr,
                                      'decay_time': decay_time,
                                      'min_SNR': min_SNR,
                                      'rval_thr': rval_thr,
                                      'use_cnn': use_cnn})

cnm.estimates.evaluate_components(images, cnm.params, dview=dview)

print(('Keeping ' + str(len(cnm.estimates.idx_components)) +
       ' and discarding  ' + str(len(cnm.estimates.idx_components_bad))))
#%%
cnm.params.set('temporal', {'p': p})
cnm2 = cnm.refit(images)