from pathlib import Path
from tempfile import TemporaryDirectory
import matplotlib.pyplot as plt
import numpy as np
import suite2p

import json
import pprint as pp
import matplotlib as mpl

# Figure Style settings for notebook.
mpl.rcParams.update({
	'axes.spines.left': False,
	'axes.spines.bottom': False,
	'axes.spines.top': False,
	'axes.spines.right': False,
	'legend.frameon': False,
	'figure.subplot.wspace': .01,
	'figure.subplot.hspace': .01,
	'figure.figsize': (18, 13),
	'ytick.major.left': False,
})
jet = mpl.cm.get_cmap('jet')
jet.set_bad(color='k')


# %%

def visualize_registrations(output_op):
	fig, axarr = plt.subplots(1, 4)

	axarr[0].imshow(output_op['refImg'], cmap='gray', )
	axarr[0].set_title("Reference Image for Registration");

	axarr[1].imshow(output_op['max_proj'], cmap='gray')
	axarr[1].set_title("Registered Image, Max Projection");

	axarr[2].imshow(output_op['meanImg'], cmap='gray')
	axarr[2].set_title("Mean registered image")

	axarr[3].imshow(output_op['meanImgE'], cmap='gray')
	axarr[3].title("High-pass filtered Mean registered image")

	return fig, axarr


def show_cell_footprints(output_stats, output_op):
	im = suite2p.ROI.stats_dicts_to_3d_array(stats, Ly=output_op['Ly'], Lx=output_op['Lx'], label_id=True)
	im[im == 0] = np.nan

	fig, axarr = plt.subplots(1, 4)

	axarr[0].imshow(output_op['max_proj'], cmap='gray')
	axarr[0].set_title("Registered Image, Max Projection")

	axarr[1].imshow(np.nanmax(im, axis=0), cmap='jet')
	axarr[1].set_title("All ROIs Found")

	axarr[2].imshow(np.nanmax(im[~iscell], axis=0, ), cmap='jet')
	axarr[2].set_title("All Non-Cell ROIs")

	axarr[3].imshow(np.nanmax(im[iscell], axis=0), cmap='jet')
	axarr[3].set_title("All Cell ROIs")


def show_roi_traces(f, f_neu, sp, rois):
	fig = plt.figure(figsize=[20, 20])

	plt.suptitle("Flourescence and Deconvolved Traces for Different ROIs", y=0.92)
	axarr = []
	for i, roi in enumerate(rois):
		axarr.append(plt.subplot(len(rois), 1, i + 1, ))
		f = f_cells[roi]
		f_neu = f_neuropils[roi]
		sp = spks[roi]
		# Adjust spks range to match range of fluroescence traces
		fmax = np.maximum(f.max(), f_neu.max())
		fmin = np.minimum(f.min(), f_neu.min())
		frange = fmax - fmin
		sp /= sp.max()
		sp *= frange
		plt.plot(f, label="Cell Fluorescence")
		plt.plot(f_neu, label="Neuropil Fluorescence")
		plt.plot(sp + fmin, label="Deconvolved")
		plt.xticks(np.arange(0, f_cells.shape[1], f_cells.shape[1] / 10))
		plt.ylabel(f"ROI {roi}", rotation=0)
		plt.xlabel("frame")
		if i == 0:
			plt.legend(bbox_to_anchor=(0.93, 2))
	return fig, axarr


# %% LOAD METADATA FROM JSON FILES IN MOVIE FOLDER

tiff_folder = Path(r"G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\data\processed_data\2021-07-27\1\movie_001")

print("____________MOVIE METADATA____________")

xml_file = tiff_folder.joinpath("xml_meta.json")
with open(xml_file, 'r') as f:
	xml_meta = json.load(f)
print(f"\nXML_META:")
print("------------")
pp.pprint(xml_meta)

sync_file = tiff_folder.joinpath("sync_meta.json")
with open(sync_file, 'r') as f:
	sync_meta = json.load(f)
print(f"\nSYNC_META:")
print("--------------")
pp.pprint(sync_meta)

# %% Make ops and db for suite2p

ops = np.load(
	Path(r"G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\data\processed_data\2021-07-27\1\movie_001\ops.npy"),
	allow_pickle=True).item()

db = dict(data_path=[str(tiff_folder)],
          tau=10.0,
          nimg_init=500,
          block_size=[56, 56],
          batch_size=xml_meta['timepoints'] / sync_meta['n_blocks'],
          nplanes=xml_meta['z_steps'],
          ignore_flyback=[],
          diameter=5,
          anatomical_only=2.0,
          )
# %% RUN SUITE2P

output_ops = suite2p.run_s2p(ops=ops, db=db)

# %%
# output_op = output_ops[0]
current_ops = output_ops[0]
results_files = list(Path(current_ops['save_path']).iterdir())
# %% visualize registrations

# %%  Load cell footprints

stats_file = Path(output_op['save_path']).joinpath('stat.npy')
iscell = np.load(Path(output_op['save_path']).joinpath('iscell.npy'), allow_pickle=True)[:, 0].astype(bool)
stats = np.load(stats_file, allow_pickle=True)
stats.shape, iscell.shape

# %% plot cell traces
f_cells = np.load(Path(output_op['save_path']).joinpath('F.npy'))
f_neuropils = np.load(Path(output_op['save_path']).joinpath('Fneu.npy'))
spks = np.load(Path(output_op['save_path']).joinpath('spks.npy'))
f_cells.shape, f_neuropils.shape, spks.shape
