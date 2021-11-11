#!/usr/bin/env python
# coding: utf-8

# ## Motion Correction demo
# 
# This notebook demonstrates the various routines for motion correction in the CaImAn package. It demonstrates the usage of rigid and piecewise rigid motion correction on a two-photon calcium imaging dataset using the NoRMCorre algorithm [[1]](#normcorre), as well as several measures for quality assessment. This notebook should be interpreted more as a tutorial of the various methods. In practice, you can use either rigid or piecewise rigid motion correction depending on the motion of the dataset.
# 
# The dataset used in this notebook is provided by Sue Ann Koay and David Tank, Princeton University. This is a two photon calcium imaging dataset. For motion correction of one photon microendoscopic data the procedure is similar, with the difference, that the shifts are inferred on high pass spatially filtered version of the data. For more information check the demos for one photon data in the CaImAn package.
# 
# More information about the NoRMCorre algorithm can be found in the following paper:
# 
# <a name="normcorre"></a>[1] Pnevmatikakis, E.A., and Giovannucci A. (2017). NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data. Journal of Neuroscience Methods, 291:83-92 [[paper]](https://doi.org/10.1016/j.jneumeth.2017.07.031)

# In[ ]:


from builtins import zip
from builtins import str
from builtins import map
from builtins import range
from past.utils import old_div

import cv2
import matplotlib.pyplot as plt
import numpy as np
import logging
from pathlib import Path
import caiman as cm
from caiman.motion_correction import MotionCorrect, tile_and_correct, motion_correction_piecewise


try:
    cv2.setNumThreads(0)
except:
    pass

try:
    if __IPYTHON__:
        get_ipython().magic('load_ext autoreload')
        get_ipython().magic('autoreload 2')
except NameError:
    pass

logging.basicConfig(format=
                    "%(relativeCreated)12d [%(filename)s:%(funcName)20s():%(lineno)s] [%(process)d] %(message)s",
                    # filename="/tmp/caiman.log",
                    level=logging.DEBUG)



# First download the file and load it in memory to view it. Note that it is not necessary to load the file in memory in order to perform motion correction. Here we load it to inspect it. Viewing the file occurs with OpenCV and will a open a new window. **To exit click on the video and press q.**
# 
# The `download_demo` function will download the specific file for you and return the complete path to the file which will be stored in your `caiman_data` directory. If you adapt this demo for your data make sure to pass the complete path to your file(s). Remember to pass the `fname` variable as a list.

# In[ ]:

PROC_DIR = Path("/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/"
                "data/processed_data/2021-08-11/1/movie_001_ds5")
fnames = sorted(PROC_DIR.glob('*.tif'))
fnames = [str(item) for item in fnames]
print(fnames)
#%%

m_orig = cm.load_movie_chain(fnames)
T, d1, d2, d3 = m_orig.shape
print(m_orig.shape)
#downsample_ratio = .2  # motion can be perceived better when downsampling in time
#m_orig.resize(1, 1, 1, downsample_ratio).play(q_max=99.5, fr=30, magnification=2)  # play movie (press q to exit)
#%%
Cn = cm.local_correlations(m_orig, swap_dim=False)
#%%
dims = Cn.shape
d1, d2, d3 = dims
# Now set some parameters that are used for motion correction.

x, y = (int(2 * (d1 + d3)), int(2 * (d2 + d3)))
scale = 6/x
fig = plt.figure(figsize=(scale*x, scale*y))
axz = fig.add_axes([1-d1/x, 1-d2/y, d1/x, d2/y])
plt.imshow(Cn.max(2).T, cmap='gray')
plt.title('Max.proj. z')
plt.xlabel('x')
plt.ylabel('y')
axy = fig.add_axes([0, 1-d2/y, d3/x, d2/y])
plt.imshow(Cn.max(0), cmap='gray')
plt.title('Max.proj. x')
plt.xlabel('z')
plt.ylabel('y')
axx = fig.add_axes([1-d1/x, 0, d1/x, d3/y])
plt.imshow(Cn.max(1).T, cmap='gray')
plt.title('Max.proj. y')
plt.xlabel('x')
plt.ylabel('z');
plt.show()

#%%
from caiman.source_extraction.cnmf import cnmf as cnmf


# motion correction parameters
opts_dict = {'fnames': fnames,
             'fr': 3.0,
             'dims': (256, 256, 16),
            'strides': (64, 64, 4),    # start a new patch for pw-rigid motion correction every x pixels
            'overlaps': (16, 16, 2),   # overlap between pathes (size of patch strides+overlaps)
            'max_shifts': (10, 10, 2),   # maximum allowed rigid shifts (in pixels)
            'max_deviation_rigid': 5,  # maximum shifts deviation allowed for patch with respect to rigid shifts
            'pw_rigid': True,         # flag for performing non-rigid motion correction
            'is3D': True,
             }

opts = cnmf.params.CNMFParams(params_dict=opts_dict)


# In[ ]:


max_shifts = (6, 6)  # maximum allowed rigid shift in pixels (view the movie to get a sense of motion)
strides = (48, 48)  # create a new patch every x pixels for pw-rigid correction
overlaps = (24, 24)  # overlap between pathes (size of patch strides+overlaps)
num_frames_split = 100  # length in frames of each chunk of the movie (to be processed in parallel)
max_deviation_rigid = 3  # maximum deviation allowed for patch with respect to rigid shifts
pw_rigid = False  # flag for performing rigid or piecewise rigid motion correction
shifts_opencv = True  # flag for correcting motion using bicubic interpolation (otherwise FFT interpolation is used)
border_nan = 'copy'  # replicate values along the boundary (if True, fill in with NaN)

# Note that here the data presented here has been downsampled in space by a factor of 2 to reduce the file size. As a
# result the spatial resolution is coarser here (around 2 microns per pixel). If we were operating at the original
# resolution, several of the parameters above, e.g., ```max_shifts, strides, overlaps, max_deviation_rigid```,
# could have been larger by a factor of 2.

# ###### Motion correction is performed in parallel on chunks taken across times.
# 
# We first need to start a cluster. The default backend mode for parallel processing is through the multiprocessing package. To make sure that this package is viewable from everywhere before starting the notebook these commands need to be executed from the terminal (in Linux and Windows):
# ```bash
#    export MKL_NUM_THREADS=1
#    export OPENBLAS_NUM_THREADS=1 
#    ```

# In[ ]:


# %% SETUP CLUSTER : start the cluster (if a cluster already exists terminate it)
if 'dview' in locals():
    cm.stop_server(dview=dview)
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=False)

# We first need to create a motion correction object with the parameters specified above. We pass directly its input
# arguments in the constructor below. Alternatively, we can use the `params` object and construct it by passing the
# arguments of `params.motion`. See the notebook `demo_pipeline.ipynb` for an example of this usage.

#%% MOCO PARAMS


# create a motion correction object
mc = MotionCorrect(fnames,
                   dview=dview,
                   max_shifts=max_shifts,
                   strides=strides,
                   overlaps=overlaps,
                   max_deviation_rigid=max_deviation_rigid,
                   shifts_opencv=shifts_opencv,
                   nonneg_movie=True,
                   border_nan=border_nan)

# <h1> Rigid motion correction</h1>
# <p> The original file exhibits a lot of motion. In order to correct for it we are first trying a simple rigid motion correction algorithm. This has already been selected by setting the parameter `pw_rigid=False` during the construction of the `MotionCorrect` object. The algorithm first creates a template by averaging frames from the video. It then tries to match each frame to this template. In addition the template will get updated during the matching process, resulting in a single precise template that is used for subpixel registration.  </p>
# <img src="../../docs/img/rigidcorrection.png" />

# In[ ]:


get_ipython().run_cell_magic('capture', '',
                             '# correct for rigid motion correction and save the file (in memory mapped form)\nmc.motion_correct(save_movie=True)')

# The motion corrected file is automatically save as memory mapped file in the location given by `mc.mmap_file`. The rigid shifts are also save in `mc.shifts_rig`.

# In[ ]:


# load motion corrected movie
m_rig = cm.load(mc.mmap_file)
bord_px_rig = np.ceil(np.max(mc.shifts_rig)).astype(np.int)
# %% visualize templates
plt.figure(figsize=(20, 10))
plt.imshow(mc.total_template_rig, cmap='gray')

# In[ ]:


# %% inspect movie
m_rig.resize(1, 1, downsample_ratio).play(
    q_max=99.5, fr=30, magnification=2, bord_px=0 * bord_px_rig)  # press q to exit

# plot the shifts computed by rigid registration

# In[ ]:


# %% plot rigid shifts
plt.close()
plt.figure(figsize=(20, 10))
plt.plot(mc.shifts_rig)
plt.legend(['x shifts', 'y shifts'])
plt.xlabel('frames')
plt.ylabel('pixels')

# ## Piecewise rigid registration
# While rigid registration corrected for a lot of the movement, there is still non-uniform motion present in the registered file. To correct for that we can use piece-wise rigid registration **directly in the original file** by setting `mc.pw_rigid=True`. As before the registered file is saved in a memory mapped format in the location given by `mc.mmap_file`.

# In[ ]:


get_ipython().run_cell_magic('capture', '',
                             '#%% motion correct piecewise rigid\nmc.pw_rigid = True  # turn the flag to True for pw-rigid motion correction\nmc.template = mc.mmap_file  # use the template obtained before to save in computation (optional)\n\nmc.motion_correct(save_movie=True, template=mc.total_template_rig)\nm_els = cm.load(mc.fname_tot_els)\nm_els.resize(1, 1, downsample_ratio).play(\n    q_max=99.5, fr=30, magnification=2,bord_px = bord_px_rig)')

# Now concatenate all the movies (raw, rigid, and pw-rigid) for inspection

# In[ ]:


cm.concatenate([m_orig.resize(1, 1, downsample_ratio) - mc.min_mov * mc.nonneg_movie,
                m_rig.resize(1, 1, downsample_ratio), m_els.resize(
        1, 1, downsample_ratio)], axis=2).play(fr=60, q_max=99.5, magnification=2, bord_px=bord_px_rig)

# From the movie we can see that pw-rigid registration corrected for the non uniform motion of the data. This was done by estimating different displacement vectors for the different patches in the FOV. This can be visualized by plotting all the computed shifts were a dispersion in the shifts in the y direction is apparent. In this case, the shifts along the two axes are stored in `mc.x_shifts_els` and `mc.y_shifts_els`, respectively.

# In[ ]:


# %% visualize elastic shifts
plt.close()
plt.figure(figsize=(20, 10))
plt.subplot(2, 1, 1)
plt.plot(mc.x_shifts_els)
plt.ylabel('x shifts (pixels)')
plt.subplot(2, 1, 2)
plt.plot(mc.y_shifts_els)
plt.ylabel('y_shifts (pixels)')
plt.xlabel('frames')
# %% compute borders to exclude
bord_px_els = np.ceil(np.maximum(np.max(np.abs(mc.x_shifts_els)),
                                 np.max(np.abs(mc.y_shifts_els)))).astype(np.int)

# The improvement in performance can also be seen by a more crisp summary statistic image. Below we plot the correlation images for the three datasets.

# In[ ]:


plt.figure(figsize=(20, 10))
plt.subplot(1, 3, 1);
plt.imshow(m_orig.local_correlations(eight_neighbours=True, swap_dim=False))
plt.subplot(1, 3, 2);
plt.imshow(m_rig.local_correlations(eight_neighbours=True, swap_dim=False))
plt.subplot(1, 3, 3);
plt.imshow(m_els.local_correlations(eight_neighbours=True, swap_dim=False))

# In[ ]:


cm.stop_server(dview=dview)  # stop the server

# ## Quality assessment
# 
# Apart from inspection, the performance of the registration methods can be quantified using several measures. Below we compute measures such as correlation of each frame with mean, crispness of summary image, and residual optical flow for all three cases. For more info see [[1]](#normcorre). Note that computation of the residual optical flow can be computationally intensive.

# In[ ]:


get_ipython().run_cell_magic('capture', '',
                             '#% compute metrics for the results (TAKES TIME!!)\nfinal_size = np.subtract(mc.total_template_els.shape, 2 * bord_px_els) # remove pixels in the boundaries\nwinsize = 100\nswap_dim = False\nresize_fact_flow = .2    # downsample for computing ROF\n\ntmpl_orig, correlations_orig, flows_orig, norms_orig, crispness_orig = cm.motion_correction.compute_metrics_motion_correction(\n    fnames[0], final_size[0], final_size[1], swap_dim, winsize=winsize, play_flow=False, resize_fact_flow=resize_fact_flow)\n\ntmpl_rig, correlations_rig, flows_rig, norms_rig, crispness_rig = cm.motion_correction.compute_metrics_motion_correction(\n    mc.fname_tot_rig[0], final_size[0], final_size[1],\n    swap_dim, winsize=winsize, play_flow=False, resize_fact_flow=resize_fact_flow)\n\ntmpl_els, correlations_els, flows_els, norms_els, crispness_els = cm.motion_correction.compute_metrics_motion_correction(\n    mc.fname_tot_els[0], final_size[0], final_size[1],\n    swap_dim, winsize=winsize, play_flow=False, resize_fact_flow=resize_fact_flow)')

# Plot correlation with mean frame for each dataset

# In[ ]:


plt.figure(figsize=(20, 10))
plt.subplot(211);
plt.plot(correlations_orig);
plt.plot(correlations_rig);
plt.plot(correlations_els)
plt.legend(['Original', 'Rigid', 'PW-Rigid'])
plt.subplot(223);
plt.scatter(correlations_orig, correlations_rig);
plt.xlabel('Original');
plt.ylabel('Rigid');
plt.plot([0.3, 0.7], [0.3, 0.7], 'r--')
axes = plt.gca();
axes.set_xlim([0.3, 0.7]);
axes.set_ylim([0.3, 0.7]);
plt.axis('square');
plt.subplot(224);
plt.scatter(correlations_rig, correlations_els);
plt.xlabel('Rigid');
plt.ylabel('PW-Rigid');
plt.plot([0.3, 0.7], [0.3, 0.7], 'r--')
axes = plt.gca();
axes.set_xlim([0.3, 0.7]);
axes.set_ylim([0.3, 0.7]);
plt.axis('square');

# In[ ]:


# print crispness values
print('Crispness original: ' + str(int(crispness_orig)))
print('Crispness rigid: ' + str(int(crispness_rig)))
print('Crispness elastic: ' + str(int(crispness_els)))

# In[ ]:


# %% plot the results of Residual Optical Flow
fls = [mc.fname_tot_els[0][:-4] + '_metrics.npz', mc.fname_tot_rig[0][:-4] +
       '_metrics.npz', mc.fname[0][:-4] + '_metrics.npz']

plt.figure(figsize=(20, 10))
for cnt, fl, metr in zip(range(len(fls)), fls, ['pw_rigid', 'rigid', 'raw']):
    with np.load(fl) as ld:
        print(ld.keys())
        print(fl)
        print(str(np.mean(ld['norms'])) + '+/-' + str(np.std(ld['norms'])) +
              ' ; ' + str(ld['smoothness']) + ' ; ' + str(ld['smoothness_corr']))

        plt.subplot(len(fls), 3, 1 + 3 * cnt)
        plt.ylabel(metr)
        try:
            mean_img = np.mean(
                cm.load(fl[:-12] + 'mmap'), 0)[12:-12, 12:-12]
        except:
            try:
                mean_img = np.mean(
                    cm.load(fl[:-12] + '.tif'), 0)[12:-12, 12:-12]
            except:
                mean_img = np.mean(
                    cm.load(fl[:-12] + 'hdf5'), 0)[12:-12, 12:-12]

        lq, hq = np.nanpercentile(mean_img, [.5, 99.5])
        plt.imshow(mean_img, vmin=lq, vmax=hq)
        plt.title('Mean')
        plt.subplot(len(fls), 3, 3 * cnt + 2)
        plt.imshow(ld['img_corr'], vmin=0, vmax=.35)
        plt.title('Corr image')
        plt.subplot(len(fls), 3, 3 * cnt + 3)
        # plt.plot(ld['norms'])
        # plt.xlabel('frame')
        # plt.ylabel('norm opt flow')
        # plt.subplot(len(fls), 3, 3 * cnt + 3)
        flows = ld['flows']
        plt.imshow(np.mean(
            np.sqrt(flows[:, :, :, 0] ** 2 + flows[:, :, :, 1] ** 2), 0), vmin=0, vmax=0.3)
        plt.colorbar()
        plt.title('Mean optical flow')
