from pathlib import PureWindowsPath, PurePath, Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from scipy import spatial, signal, stats, ndimage, optimize

from collections import namedtuple
from functools import partial
import copy
import json
import itertools

from sklearn import preprocessing, svm, model_selection, decomposition, pipeline

# from suite2p.extraction import dcnv
import caiman
import caiman.source_extraction.cnmf as cm

from ryim import dirutils
from rearing.ca import respvec, vis
import rearing_datasets
from step00_process_suite2p import load_metadata, get_title_str, postprocess_metadata

import importlib as imp

mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 10,
                     'savefig.format': 'png',
                     'figure.titlesize': 'medium',
                     'axes.labelsize': 'small',
                     'axes.titlesize': 'small',
                     'xtick.labelsize': 'x-small',
                     'ytick.labelsize': 'x-small',
                     'legend.fontsize': 'x-small'})


# %%

def select_cells(cellprob, iplane, cell_thr, select_plane):
	cellthr_idx = (cellprob > cell_thr)
	plane_idx = np.isin(iplane, select_plane)
	good_cell_idx = cellthr_idx & plane_idx
	return good_cell_idx


def get_margin(coef):
	margin = 1 / np.sqrt(np.sum(coef ** 2))
	# margin = 1 / np.sum(np.abs(coef))
	return margin


def make_ovo_training_data(trials_x_cells, df_stimulus, odor_a, odor_b):
	"""
	Converts peak_amp, df_stimulus into format that sklearn models can use.

	ovo = one vs one

	:param trials_x_cells: peak measurement array, transposed (each trial is a row, each cell a column)
	:type trials_x_cells: np.ndarray
	:param df_stimulus:
	:param odor_a:
	:param odor_b:
	:return:
		- X: samples x features [# of trials x n_cells]
		- y: labels (odor a = 0, odor_b = 1)
		- le: label encoder (sklearn.preprocessing)
	"""

	mask = df_stimulus['stimname'].isin([odor_a, odor_b])
	X = trials_x_cells[mask, :]
	y = df_stimulus.loc[mask, 'stimname'].to_numpy()
	le = preprocessing.LabelEncoder().fit([odor_a, odor_b])
	return X, y, le


def make_AB_detector_training_data(trials_x_cells, df_stimulus, odor_ab):
	"""
	:param trials_x_cells: peak measurement array, transposed (each trial is a row, each cell a column)
	:type trials_x_cells: np.ndarray
	:param df_stimulus:
	:param odor_ab:
	:return:
		- X: samples x features [# of trials x n_cells]
		- y: labels (odor a = 0, odor_b = 1)

	Converts peak_amp, df_stimulus into format that sklearn models can use.
	"""

	mask = df_stimulus['stimname'] == odor_ab
	X = trials_x_cells[mask, :]
	y = df_stimulus.loc[mask, 'stimname'].to_numpy()
	return X, y


def tidy_proj_mat(splits_x_trials, df_stimulus):
	"""
	Turns Z_mat into long form/tidy for easy plotting
	:param splits_x_trials:
	:param df_stimulus:
	:return:
	"""
	df = pd.DataFrame(splits_x_trials).melt().set_index('variable')
	df = df.join(df_stimulus)
	return df


def gather_svm_results(dataset_list, svm_filename):
	"""
	Collects svm metadata and svm results into lists
		- loads results of one set of parameters (defined by svm_filename)
		- loads data for multiple flies & movies into a flat list

	:param dataset_list: List of namedtuple Datasets (see rearing_datasets.py)
	:param svm_filename: which file in suite2p/**/step03_svm_analysis to load (based off parameters)
						   ex - norml1__c10__svc
	:return: dict w/ keys
				ainfo_list,
				title_str_list,
				svm_analysis_list,
				tidy_Z_list,
				param_dict_list
	"""
	ainfo_list = []
	title_str_list = []
	svm_analysis_list = []
	tidy_Z_list = []
	param_dict_list = []

	for i, ds in enumerate(dataset_list):
		stat_files = rearing_datasets.path_to_statfiles(ds)

		for f in stat_files:
			ainfo = dirutils.parse_dataset_info(f)

			nt_meta_data, title_str, peak_data, oasis_results, iscell = extract(f)

			loadfile = f.parent.joinpath('step03_svm_analysis', f"{svm_filename}.npy")
			svm_data = np.load(loadfile, allow_pickle=True).item()

			# add columns to tidy_Z dataframe
			tidy_Z = svm_data['tidy_Z']
			tidy_Z['ifly'] = i
			tidy_Z['flydate'] = ainfo.flydate
			tidy_Z['flynum'] = ainfo.flynum
			tidy_Z['rear_type'] = nt_meta_data.rear_type
			tidy_Z['blk_conc'] = nt_meta_data.blk_conc

			analysis_dict = dict(ainfo=ainfo,
			                     nt_meta_data=nt_meta_data,
			                     title_str=title_str,
			                     peak_data=peak_data,
			                     iscell=iscell,
			                     oasis_results=oasis_results, )
			ainfo_list.append(ainfo)
			title_str_list.append(title_str)
			svm_analysis_list.append(svm_data)
			tidy_Z_list.append(tidy_Z)
			param_dict_list.append(svm_data['param_dict'])

	svm_results = dict(ainfo_list=ainfo, title_str_list=title_str_list, svm_analysis_list=svm_analysis_list,
	                   tidy_Z_list=tidy_Z_list, param_dict_list=param_dict_list)
	# return ainfo_list, title_str_list, svm_analysis_list, tidy_Z_list, param_dict_list

	return svm_results


def combine_tidy_Z_list(tidy_Z_list):
	"""
	Combines list of tidy_Z dataframes (Z_mat in long format) into a single pd.DataFrame,
	and computes means for groups defined by ['ifly', 'blk_conc', 'stimname'].

	:param tidy_Z_list: list of tidy_Z dataframes
						see 'tidy_Z_list' in dict returned by function gather_svm_results(...)

	:type tidy_Z_list: list[pd.DataFrame]
	:return: tidy_Z_cat, tidy_Zmean_cat in dict
	:rtype: dict

	"""
	tidy_Z_cat = pd.concat(tidy_Z_list, axis=0)
	tidy_Zmean_cat = tidy_Z_cat.groupby(by=['ifly', 'blk_conc', 'stimname']).mean()
	tidy_Zmean_cat.reset_index(inplace=True)

	return dict(tidy_Z_cat=tidy_Z_cat, tidy_Zmean_cat=tidy_Zmean_cat)


def extract(f):
	"""
	Load nt_meta_data, title string, and peak amplitude matrix from dataset

	:param f: filepath to stat.npy
	:type f: Path | str
	:return:
		-   nt_meta_data: namedtuple with the fields
							df_stimulus",
						   "df_frametimeinfo",
						   "timestamps",
						   "session_info",
						   "xml_meta",
						   "sync_meta",
						   "meta_yaml",
						   "fs", "blk_conc", "rear_type"

	   -      title_str: string about dataset, for plot titles
						(date, flynum, movie, block, concentration, rearing condition)

	   - cells_x_trials: dict of arrays w/ 1 column per trial/stimulus presentation and 1 row per cell,
								and respvec_ops (parameters for peak quantification)
						 contains some peak quantification value, loaded from suite2p/**/step01_extract_peaks

	:rtype: namedtuple. str, dict
	"""
	ainfo = dirutils.parse_dataset_info(f)
	nt_meta_data = postprocess_metadata(load_metadata(f))
	ainfo = dirutils.parse_dataset_info(f)
	title_str = get_title_str(ainfo, nt_meta_data)

	folder = f.parent.joinpath('step01_extract_peaks')
	baseline_mean = np.load(folder.joinpath('baseline_mean.npy'))
	baseline_std = np.load(folder.joinpath('baseline_std.npy'))
	peak_amp = np.load(folder.joinpath('peak_amp.npy'))
	peak_mean = np.load(folder.joinpath('peak_mean.npy'))
	respvec_ops = np.load(folder.joinpath('respvec_ops.npy'), allow_pickle=True)

	peak_data = dict(baseline_mean=baseline_mean,
	                 baseline_std=baseline_std,
	                 peak_amp=peak_amp,
	                 peak_mean=peak_mean,
	                 respvec_ops=respvec_ops)

	folder = f.parent.joinpath('step00_process_suite2p')
	oasis_results = np.load(folder.joinpath('oasis_results.npy'), allow_pickle=True)

	iscell = np.load(f.with_name('iscell.npy'))
	n_cells = iscell.shape[0]

	if f.parent.name == 'combined':
		iplane = np.load(f.with_name('iplane.npy')).squeeze()
		iscell = np.append(iscell, iplane[:, np.newaxis], axis=1)
	else:
		iscell = np.append(iscell, np.zeros((n_cells, 1)), axis=1)

	return nt_meta_data, title_str, peak_data, oasis_results, iscell


def transform(X0, df_stimulus, odor_a, odor_b, param_dict, verbose=True):
	"""
	Runs SVM (odor_a vs odor_b) classifier pipeline for multiple splits of X0, w/ classifier determined by param_dict.

	Returns list of trained classifiers, classifier performance, and the decision function output for X0.

	:param X0:
	:type X0:
	:param df_stimulus:
	:type df_stimulus:
	:param odor_a:
	:type odor_a:
	:param odor_b:
	:type odor_b:
	:param param_dict:
	:type param_dict:
	:return:
	:rtype: dict w/ keys
		-               X0 : training data, generated from peaks
		-       param_dict :  list of parameters to run through
								- norm: 'l1' or 'l2', input to preprocessing.Normalizer
								- svm_type: 'svm' or 'linear_svc'
								- c: regularization parameter
		-              ss : modelSelection.ShuffleSplit
		- train_index_list: list of lists
							train_index values returned by each iteration of ss.split(X)
		- train_index_list: list of lists,
							train_index values returned by each iteration of ss.split(X)
		-         clf_list: list of fit classifiers
		-      score_list : list of performance scores for each classifier in clf_list

								ex:  [ clf.score(X_test, y_test) for clf in clf_list ]

	   -            tidy_Z: decision function values for Z_mat transformed into long-form,
							where Z_mat is a [splits x trial] np.array of clf.decision_function values

								ex: Z_mat = [clf.decision_function[X0] for clf in clf_list]
									tidy_Z = tidy_proj_mat(Z_mat)
		-          verbose: whether or not to print info
	"""

	print(param_dict)

	X, y, le = make_ovo_training_data(X0, df_stimulus, odor_a, odor_b)

	train_index_list = []
	test_index_list = []
	clf_list = []
	score_list = []
	Z_list = []
	coef_list = []

	ss = model_selection.ShuffleSplit(n_splits=10, test_size=0.3, random_state=0)

	for train_index, test_index in ss.split(X):
		X_train = X[train_index, :]
		y_train = y[train_index]
		X_test = X[test_index, :]
		y_test = y[test_index]
		train_index_list.append(train_index)
		test_index_list.append(test_index)

		if verbose:
			print("%s %s" % (train_index, test_index))

		if param_dict['svm_type'] == 'svc':
			clf = pipeline.make_pipeline(preprocessing.Normalizer(norm=param_dict['norm']),
			                             svm.SVC(kernel='linear',
			                                     class_weight='balanced',
			                                     C=param_dict['c'],
			                                     probability=True,
			                                     verbose=verbose,
			                                     decision_function_shape='ovo')
			                             )

		elif param_dict['svm_type'] == 'linear_svc':
			clf = pipeline.make_pipeline(preprocessing.Normalizer(norm=param_dict['norm']),
			                             svm.LinearSVC(penalty='l1',
			                                           class_weight='balanced',
			                                           C=param_dict['c'], )
			                             )

		clf.fit(X_train, y_train)
		clf_list.append(clf)

		score = clf.score(X_test, y_test)
		score_list.append(score)

		if verbose:
			print(f'score: {score}')

		Z = clf.decision_function(X0)
		Z_list.append(Z)

		coef_list.append(clf[param_dict['svm_type']].coef_.squeeze())

	# converts decision functions to long form
	Z_mat = np.array(Z_list)
	tidy_Z = tidy_proj_mat(Z_mat, df_stimulus)

	svm_results = dict(X0=X0, param_dict=param_dict, ss=ss,
	                   train_index_list=train_index_list,
	                   test_index_list=test_index_list,
	                   clf_list=clf_list, score_list=score_list,
	                   tidy_Z=tidy_Z)
	return svm_results


def save(outputdir, svm_results):
	"""
	Saves svm results for one set of parameters.

	:param outputdir: folder inwhich to save the .npy files
	:param svm_results:
	:return:
	"""
	param_dict = svm_results['param_dict']
	savename = f"norm{param_dict['norm']}__c{param_dict['c']}__{param_dict['svm_type']}.npy"
	np.save(outputdir.joinpath(savename), param_dict, allow_pickle=True)
	return True


def plot_repeated_decision_functions(df, ax=None):
	"""
	Plot
	:param df: tidy_Z
	:param ax:
	:return:
	"""
	if ax is None:
		fig, ax = plt.subplots()

	sns.boxplot(data=df, x='stimname', y='value',
	            width=0.5,
	            hue='stimname',
	            ax=ax,
	            palette=rearing_datasets.clr_dict)

	sns.stripplot(data=df, x='stimname', y='value', hue='stimname', ax=ax,
	              jitter=0.05,
	              palette=rearing_datasets.clr_dict, size=3)

	for patch in ax.artists:
		r, g, b, a = patch.get_facecolor()
		patch.set_facecolor((r, g, b, .3))
	plt.subplots_adjust(right=0.75, hspace=0.25)
	sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), title=None, frameon=False)

	ax.set_title('Linear classifier')
	ax.set_ylim([-1.3, 1.3])
	return fig, ax


def main(f, cell_thr_type='quantile', cell_thr=0.5, select_plane=None):
	"""
	STEP 01: LOAD DATA

	:param select_plane:
	:type select_plane:
	:param cell_thr:
	:type cell_thr:
	:param cell_thr_type:
	:type cell_thr_type:
	:param f: pathlib.Path, path to suite2p/**/stat.npy file

	:return:
	"""
	# get metadata, load suite2p outputs
	# ainfo = AnalysisSet(flydate, flynum, mov, is_blk, blk)
	plt.close('all')

	ainfo = dirutils.parse_dataset_info(f)

	# extract
	# ------------------------------------------
	nt_meta_data, title_str, peak_data, oasis_results, iscell = extract(f)

	iplane = iscell[:, 2].astype('int')
	cellprob = iscell[:, 1]
	iscell = iscell[:, 0]

	cell_thr = np.quantile(cellprob, .5).squeeze()
	select_plane = int(stats.mode(iplane).mode[0])

	good_cell_idx = select_cells(cellprob, iplane, cell_thr, select_plane)
	good_cells = np.squeeze(np.argwhere(good_cell_idx))

	peak_amp = copy.deepcopy(peak_data['peak_amp'])
	X0 = peak_amp[good_cell_idx, :].T

	# parameters to try
	grid = dict(svm_type=['svc'],
	            norm=['l1', 'l2'],
	            c=[1, 3, 5, 10])
	param_list = list(model_selection.ParameterGrid(grid))

	outputdir = f.parent.joinpath('step03_svm_analysis')
	outputdir.mkdir(parents=True, exist_ok=True)

	for param_dict in param_list:
		# transform
		print(param_dict)
		svm_data = transform(X0, nt_meta_data.df_stimulus, '1-6ol', 'EP', param_dict)
		svm_data['good_cell_idx'] = good_cell_idx
		print(svm_data.keys())

		savename = f"norm{param_dict['norm']}__c{param_dict['c']}__{param_dict['svm_type']}.npy"
		np.save(outputdir.joinpath(savename), svm_data, allow_pickle=True)
		# save
		# ------------------------------------------
		# for key, value in svm_data.items():
		#   np.save(outputdir, svm_data)
		print(f'SVM results saved to {outputdir}')
	return True


# %%

if __name__ == '__main__':
	pfo_datasets = rearing_datasets.get_pfo_datasets()
	odor_datasets = rearing_datasets.get_odor_datasets()

	pfo_stat_files = rearing_datasets.get_pfo_statfiles()
	odor_stat_files = rearing_datasets.get_odor_statfiles()

	run_svm_analysis = True
	if run_svm_analysis:
		for f in pfo_stat_files:
			main(f)

		for f in odor_stat_files:
			main(f)

# %% CHOOSE WHICH SET OF PARAMETERS TO PLOT RESULTS FOR

FIG_DIR = Path.cwd().joinpath('doc', 'figures', 'svm analysis')

svm_filename = 'norml2__c1__svc'

# %% SECTION 1.0: LOAD AND FORMAT SVM ANALYSIS RESULTS

# paraffin-reared flies
pfo_reared_results = gather_svm_results(pfo_datasets, svm_filename)
pfo_reared_results.update(combine_tidy_Z_list(pfo_reared_results['tidy_Z_list']))

# odor-reared flies
odor_reared_results = gather_svm_results(odor_datasets, svm_filename)
odor_reared_results.update(combine_tidy_Z_list(odor_reared_results['tidy_Z_list']))

# %% SECTION 2: MAKE 1 PLOT OF SVM DECISION FUNCTIONS FOR EACH DATASET
#       - each movie/chunk, should only have odors @ 1 concentration)
#       - save all plots to a single pdf (1 pdf per rearing condition)
#       - save pdf to mb_odor_rearing/doc/figures/svm analysis/<rear_type>_svm_plots_<svm_filename>

with PdfPages(FIG_DIR.joinpath(f'pfo_svm_plots_{svm_filename}.pdf')) as pdf:
	for i, (df, param_dict, title_str) in enumerate(zip(pfo_reared_results['tidy_Z_list'],
	                                                    pfo_reared_results['param_dict_list'],
	                                                    pfo_reared_results['title_str_list'])):
		print(i)
		fig, ax = plot_repeated_decision_functions(df)
		ax.set_title(title_str + f"""\n{str(param_dict)}""")
		pdf.savefig(fig)

# save all ODOR-REARED svm results to a single pdf
with PdfPages(FIG_DIR.joinpath(f'odor_svm_plots_{svm_filename}.pdf')) as pdf:
	for i, (df, param_dict, title_str) in enumerate(zip(odor_reared_results['tidy_Z_list'],
	                                                    odor_reared_results['param_dict_list'],
	                                                    odor_reared_results['title_str_list'])):
		print(i)
		fig, ax = plot_repeated_decision_functions(df)
		ax.set_title(title_str + f"""\n{str(param_dict)}""")
		pdf.savefig(fig)

# %% SECTION 3: SVM DATA PER FLY
#      - one subplot per fly, 1 figure each for pfo and odor-reared groups
#      - each subplot contains list of concentrations (-5, -4, -3)
#      - each subplot shows data from all movies for that fly
#      - save pdf to mb_odor_rearing/doc/figures/svm analysis/<rear_type>_svm_per_fly_<svm_filename>
#
#
#            odor-reared                          pfo-reared
#     + ----- + ----- + ----- +           + ----- + ----- + ----- +
#     | fly 1 | fly 2 | fly 3 |           | fly 1 | fly 2 | fly 3 |
#     + ----- + ----- + ----- +           + ----- + ----- + ----- +
#     | fly 4 | ...                       | fly 4 | ...
#     + ----- +                           + ----- +
#
#    individual subplot contains:
#                   _____________
#                  |             |
#                  | fly 1, pfo  |
#                  |             |
#                  |_____________|
#                   |||
#        odors:      A    B    C  ...
#                    ├── -3
#        conc:       ├── -4
#                    └── -5
#

# paraffin-reared
# ---------------
g_pfo = sns.catplot(data=pfo_reared_results['tidy_Z_cat'],
                    x='stimname',
                    y='value',
                    hue='blk_conc',
                    col='ifly', col_wrap=2,
                    kind='strip', jitter=0.2,
                    dodge=True,
                    linewidth=0.5,
                    height=4, aspect=1.0,
                    legend_out=True,
                    s=3, alpha=0.75)
g_pfo.set(ylim=[-1.5, 1.5])
plt.subplots_adjust(top=0.85)
g_pfo.fig.suptitle('paraffin-reared flies: svm decision function for classifier trained on 1-6ol vs EP trials,'
                   f'\n{svm_filename}'
                   '\n10 splits per movie')
plt.show()
g_pfo.fig.savefig(FIG_DIR.joinpath(f"pfo_svm_per_fly_{svm_filename}.png"))

# odor-reared
# -----------
g_odor = sns.catplot(data=odor_reared_results['tidy_Z_cat'],
                     x='stimname',
                     y='value',
                     hue='blk_conc',
                     col='ifly', col_wrap=2,
                     kind='strip', jitter=0.2,
                     dodge=True,
                     linewidth=0.75,
                     height=4, aspect=1,
                     legend_out=True,
                     s=3, alpha=0.75)
g_odor.set(ylim=[-1.5, 1.5])
plt.subplots_adjust(top=0.85)
g_odor.fig.suptitle('odor-reared flies: svm decision function for classifier trained on 1-6ol vs EP trials,'
                    f'\n{svm_filename}'
                    '\n10 splits per movie')
plt.show()
g_odor.fig.savefig(FIG_DIR.joinpath(f"odor_svm_per_fly_{svm_filename}.png"))
# g = sns.catplot(data=tidy_Z_pfo, x='stimname', y='value', hue='blk_conc', col='ifly', kind='strip', dodge=True,
#                 height=4, aspect=0.8, legend_out=True, col_wrap=3,

# %% SECTION 4: SUMMARY PLOTS OF MEAN DECISION FUNCTION VALUES
#       - mean fly/odor/concentration value calculated
#       - single subplot axis, with all concentrations shown
#       - 1 figure per rearing condition
#       - saved to [ pfo_svm_summary_{svm_filename}.pdf ] or [ odor_svm_summary_{svm_filename}.pdf ]
#               {svm_filename}_svm_summary_by_fly?
#
#                    pfo_reared              odor_reared
#                   _____________            _____________
#                  | *         + |          | *           |        * = -3
#                  | x      +    |          | x   ...     |        x = -4
#                  | +   +       |          | +           |        + = -5
#                  |_____________|          |_____________|
#                     a  b  c  d               a  b  c  d
#                       ^ odors                  ^ odors
plt.close('all')

# paraffin-reared
# ---------------
fig_pfo_Zmean, ax = plt.subplots()
sns.pointplot(data=pfo_reared_results['tidy_Zmean_cat'],
              x="stimname",
              y="value",
              hue='blk_conc',
              dodge=.25, scale=0.5, ax=ax)
plt.subplots_adjust(top=0.85)
ax.set_title('Paraffin-reared flies: svm decision function for classifier trained on 1-6ol vs EP trials,'
             '\nEach point = mean fly/odor/concentration value'
             f'\n{svm_filename}'
             '\n10 splits per movie')
savename = FIG_DIR.joinpath(f"pfo_svm_summary_{svm_filename}.pdf")
fig_pfo_Zmean.savefig(savename)

# odor-reared
fig_odor_Zmean, ax = plt.subplots()
sns.pointplot(data=odor_reared_results['tidy_Zmean_cat'],
              x="stimname",
              y="value",
              hue='blk_conc',
              dodge=.25, scale=0.5, ax=ax)
plt.subplots_adjust(top=0.85)
ax.set_title('Paraffin-reared flies: svm decision function for classifier trained on 1-6ol vs EP trials,'
             '\nEach point = mean fly/odor/concentration value'
             f'\n{svm_filename}'
             '\n10 splits per movie')
savename = FIG_DIR.joinpath(f"odor_svm_summary_{svm_filename}.png")
fig_pfo_Zmean.savefig(savename)
# -----------
# %% SECTION 4: SUMMARY PLOT, BOTH REARING GROUPS POOLED
#               - shows both pfo- and odor-reared groups on same subplots
# 				- 1 hue per rearing condition
#               - 1 subplot per odor
#               - all concentrations shown in each subplot
#                   _____________                               +++++++++++++++++++++++++++++++++++++++
#    [* = pfo ]    |             |                              +                                     +
#    [x = odor]    | odor vs pfo |  ... x [ # of odors ] --->   +   {svm_filename}_svm_summary_all    +
#                  |             |                              +++++++++++++++++++++++++++++++++++++++
#                  |_____________|
#               [odorA @ -3, -4, -5]

pfo_reared_results['tidy_Z_cat']['rear_type'] = 'pfo'
odor_reared_results['tidy_Z_cat']['rear_type'] = 'odor'

# tidy_Z_pfo['rear_type'] = 'pfo'
# tidy_Z_odor['rear_type'] = 'odor'
tidy_Z_all = pd.concat([pfo_reared_results['tidy_Z_cat'],
                        odor_reared_results['tidy_Z_cat']],
                       axis=0)

g_all_boxen = sns.catplot(data=tidy_Z_all,
                          x='blk_conc',
                          y='value',
                          hue='rear_type',
                          col='stimname', col_wrap=3,
                          kind='boxen',
                          dodge=True,
                          height=4, aspect=0.8,
                          legend_out=True)
plt.subplots_adjust(top=0.85)
g_all_boxen.fig.suptitle('SVM decision function for classifier trained on 1-6ol vs EP trials,'
                         '\n 1 point = result of 1 classifier split'
                         '\n10 splits per movie'
                         f'\n{svm_filename}')
plt.show()

savename = FIG_DIR.joinpath(f"pfo_vs_odor_reared_svm_summary_{svm_filename}.png")
g_all_boxen.fig.savefig(savename)

#%%
flierprops = dict(marker='o', markerfacecolor='None', markersize=10,  markeredgecolor='black')
# sns.boxplot(y=df.Column,orient="v",flierprops=flierprops)

g = sns.catplot(data=pfo_reared_results['tidy_Z_cat'],
                x='blk_conc',
                y='value',
                hue='ifly',
                col='stimname', col_wrap=3,
                kind='boxen',
                # scale='linear',
                linewidth=0.5,
                dodge=True,
                height=6, aspect=1.0,
                width=0.3,
                legend_out=True,)
plt.subplots_adjust(top=0.85)
g.fig.suptitle(f'pfo-reared: \n{svm_filename}')
plt.show()
#%%
#%%
fig, ax = plt.subplots()
sns.stripplot(data=tidy_Z_all,
              x='stimname',
              y='value',
              hue='ifly', ax=ax, size=3)

plt.show()
# %%
svm_filename = 'norml2__c1__svc'

tidy_Z_odor = []
ainfo_list = []
title_str_list = []
param_dict_list = []

for i, ds in enumerate(odor_datasets):
	stat_files = rearing_datasets.path_to_statfiles(ds)

	for f in stat_files:
		ainfo = dirutils.parse_dataset_info(f)
		ainfo_list.append(ainfo)

		nt_meta_data, title_str, peak_data, oasis_results, iscell = extract(f)

		loadfile = f.parent.joinpath('step03_svm_analysis', f"{svm_filename}.npy")
		svm_data = np.load(loadfile, allow_pickle=True).item()
		tidy_Z = svm_data['tidy_Z']
		tidy_Z['ifly'] = i
		tidy_Z['rear_type'] = nt_meta_data.rear_type
		tidy_Z['blk_conc'] = nt_meta_data.blk_conc

		tidy_Z_odor.append(tidy_Z)
		title_str_list.append(title_str)
		param_dict_list.append(svm_data['param_dict'])
# %%
with PdfPages(f'odor_svm_plots_{svm_filename}.pdf') as pdf:
	for i, df in enumerate(tidy_Z_odor):
		print(i)
		fig, ax = plot_repeated_decision_functions(df)
		ax.set_title(title_str_list[i] + f"""\n{str(param_dict_list[i])}""")
		pdf.savefig(fig)
# %%

# tidy_Z_odor = pd.concat(tidy_Z_odor, axis=0)

g = sns.catplot(data=pfo_reared_results['tidy_Z_cat'],
                x='stimname',
                y='value',
                hue='blk_conc',
                col='ifly',
                kind='boxen',
                dodge=True,
                legend_out=True,
                height=4, aspect=0.8, col_wrap=3,
                width=0.5)
g.set(ylim=[-1.5, 1.5])
plt.subplots_adjust(top=0.85)
g.fig.suptitle('Odor-reared flies: svm decision function for classifier trained on 1-6ol vs EP trials,'
               f'\n{svm_filename}'
               '\n10 splits per movie')
plt.show()

# %%
plt.close('all')
tidy_Zmean_odor = tidy_Z_odor.groupby(by=['ifly', 'blk_conc', 'stimname']).mean()
tidy_Zmean_odor.reset_index(inplace=True)
fig, ax = plt.subplots()
sns.pointplot(data=tidy_Zmean_cat, x="stimname", y="value", hue='blk_conc', dodge=.25, scale=0.5, ax=ax)
plt.subplots_adjust(top=0.85)
ax.set_title('Odor-reared flies: svm decision function for classifier trained on 1-6ol vs EP trials,'
             '\nEach point = mean fly/odor/concentration value'
             # f'\n{str(param_dict_list[0])}'
             f'\n{svm_filename}'
             '\n10 splits per movie')
plt.show()

# %%


for patch in ax.artists:
	r, g, b, a = patch.get_facecolor()
	patch.set_facecolor((r, g, b, .3))
plt.subplots_adjust(right=0.75, hspace=0.25)
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), title=None, frameon=False)

ax.set_title('Linear classifier')
ax.set_ylim([-1.3, 1.3])

# %%
# %%
tidy_Z_pfo['rear_type'] = 'pfo'
tidy_Z_odor['rear_type'] = 'odor'
tidy_Z_cat = pd.concat([tidy_Z_pfo, tidy_Z_odor], axis=0)
# %%
g = sns.catplot(data=tidy_Z_pfo, x='blk_conc', y='value', hue='ifly', col='stimname', kind='boxen', dodge=True,
                height=4, aspect=0.8, legend_out=True, col_wrap=3)
plt.subplots_adjust(top=0.85)
g.fig.suptitle(f'pfo-reared: \n{str(param_dict_list[0])}')
# %%
g = sns.catplot(data=tidy_Z_odor, x='blk_conc', y='value', hue='ifly', col='stimname', kind='boxen', dodge=True,
                height=4, aspect=0.8, legend_out=True, col_wrap=3)
plt.subplots_adjust(top=0.85)
g.fig.suptitle(f'odor-reared: \n{str(param_dict_list[0])}')
# %%

# %%
# #%%
#     # dist_all, df_peakresp = respvec.compute_corrmats(peak_amp, nt_meta_data.df_stimulus)
#     # fig_corrmat_all, axarr, haxarr = vis.corrmats(dist_all)
#     # fig_corrmat_all.suptitle(title_str + f"""\n# cells={X0.shape[1]}""")
#
#     # %% SVM classifier
#     # X_train, X_test, y_train, y_test = train_test_split(X, le.transform(y), test_size=0.33)
#
#
#
#     # grid = [{'norm': ['l1', 'l2'],
#     #         'c': [1, 2, 3, 5, 10]}]
#
#
#     svm_type = 'svc'
#
#
#
#     # plot decision function magnitudes
#     # ==========================================================
#     tidy_Z = tidy_proj_mat(Z_mat, nt_meta_data.df_stimulus)
#     fig, ax = plot_repeated_decision_functions(tidy_Z, nt_meta_data.df_stimulus)
#     ax.set_title(str(param_dict))
#     # sns.heatmap(np.array(coef_list), square=False, center=0, cmap='RdBu_r')
#     fig.suptitle(title_str)
#     plt.show()
#
#
#
#     #%% plot all repetitions of the ShuffleSplits
#     # ==========================================================
#     g = sns.catplot(data=df, x='reps', y='value', hue='stimname', col='stimname', kind='strip', dodge=False,
#                     height=5, aspect=0.2, legend_out=True)
#     # fig, ax = plt.subplots()
#     # sns.swarmplot(data=df, x='stimname', y='value', hue='stimname')
#     g.set_titles("{col_name}")
#     g.set_axis_labels("rep", "decision fcn")
#     g.fig.suptitle(
#         title_str + """\nDecision function for peak amp. of deconvolved traces\n1-hexanol vs. ethyl propionate""",
#         fontsize=10)
#     plt.subplots_adjust(wspace=0.2, top=0.75)
#     sns.despine(trim=True)
#
#
#
#
#     # Plot SVM decision function results
#     # ==========================================================
#     df_svm = copy.deepcopy(nt_meta_data.df_stimulus)
#     df_svm['Z'] = Z
#
#     g = sns.catplot(data=df_svm, x='reps', y='Z', hue='stimname', col='stimname', kind='bar', dodge=False,
#                     height=5, aspect=0.2, legend_out=True)
#     g.set_titles("{col_name}")
#     g.set_axis_labels("rep", "decision fcn")
#     g.fig.suptitle(
#         title_str + """\nDecision function for peak amp. of deconvolved traces\n1-hexanol vs. ethyl propionate""",
#         fontsize=10)
#     plt.subplots_adjust(wspace=0.2, top=0.75)
#     sns.despine(trim=True)
#
#     # change opacity
#     for ax in g.axes.flat:
#         for barcontainer in ax.containers:
#             for bar in barcontainer:
#                 bar.set_alpha(0.5)
#     plt.show()
#
#     #%%
#     # X, y, le, clf, score,
#     # %%
#     """
#     STEP 02: process calcium responses, make them into population vectors
#     """
#
#     # %% extract loaded files
#     df_stimulus = copy.deepcopy(nt_meta_data.df_stimulus)
#     df_frametimeinfo = copy.deepcopy(nt_meta_data.df_frametimeinfo)
#     timestamps = copy.deepcopy(nt_meta_data.timestamps)
#     ti = copy.deepcopy(nt_meta_data.sync_meta)
#     iscell = nt_s2p_data.iscell[:, 0].astype(np.bool_)
#     cellprob = nt_s2p_data.iscell[:, 1]
#
#     if nt_meta_data.xml_meta['z_steps'] == 1:
#         use_saved_cells = False
#     else:
#         use_saved_cells = False
#
#     if use_saved_cells:
#         good_cells_bool = iscell
#     else:
#         good_cells_bool = (cellprob > np.quantile(cellprob, .5)).squeeze()
#     good_cells = np.argwhere(good_cells_bool).squeeze()
#
#     # %% compute quantile-corrected cell traces, cell ranking by probability
#     Fc = nt_s2p_data.F - 0.7 * nt_s2p_data.Fneu
#     n_cells, T = Fc.shape
#     cell_ranking = np.argsort(cellprob)[::-1]
#
#     # baseline correction w/ suite2p module
#     if nt_meta_data.xml_meta['z_steps'] > 1:
#         Fcc = dcnv.preprocess(F=Fc, baseline='maximin', win_baseline=45, sig_baseline=3, fs=nt_meta_data.fs,
#                               prctile_baseline=5)
#     else:
#         Fc = trace.split_traces_to_blocks(Fc, nt_meta_data.df_frametimeinfo.scope_pulse_idx)
#         Fcc = [dcnv.preprocess(F=item, baseline='maximin', win_baseline=60, sig_baseline=3, fs=nt_meta_data.fs,
#                                prctile_baseline=20) for item in Fc]
#         Fcc = np.hstack(tuple(Fcc))
#
#     # quantile matching & transform, quantile_range=(0.25, 0.75)
#     d = preprocessing.RobustScaler().fit_transform(Fcc)
#
#     # %% multipage pdf of all the fluorescence traces,
#     save_stacked_traces = False
#     if save_stacked_traces:
#         with PdfPages('stacked_traces.pdf') as pdf:
#             for cid0 in np.arange(0, int(np.ceil(good_cells.size / 20) * 20), 20):
#                 cid_idx = np.arange(20) + cid0
#                 cids = list(good_cells[cid_idx])
#                 fig, axarr = vis.blocked_stackplot(d, timestamps, cids=cids,
#                                                    df_stimulus=df_stimulus,
#                                                    df_frametimeinfo=df_frametimeinfo,
#                                                    peak_start=0.5, peak_len=5)
#                 plt.show(block=False)
#                 plt.pause(0.01)
#                 pdf.savefig(figure=fig)
#                 plt.close(fig)
#
#
#     # %% plot oasis results
#     save_oasis_plots = False
#
#     if save_oasis_plots:
#         # of cells per figure
#         n_axes = 40
#         for cid0 in range(0, n_axes, 41):
#             fig, axarr = plt.subplots(10, int(np.ceil(n_axes / 10)), figsize=(11, 8.5),
#                                       sharey='all', sharex='all')
#
#             for cid, ax in zip(range(cid0, cid0 + n_axes), axarr.flatten()):
#                 ax = draw_oasis_results(d[cid, :], C_dec[cid, :], Sp[cid, :],
#                                         sn_=Sn[cid], x_=timestamps, ax_=ax)
#
#                 ax.set_title(f"cell # {cid} (snr={round(Sn[cid], 2)})", fontsize=8)
#                 ax.spines[["left", "bottom"]].set_position(("data", 0))
#                 ax.set_ylim([-2, 10])
#
#         plt.subplots_adjust(bottom=0.05, top=0.9, hspace=0.25)
#         plt.suptitle(title_str, fontsize=10)
#         plt.show(block=False)
#
#     # %% compute response magnitudes
#     respvec_ops = dict(peak_start=0.5, peak_len=5, baseline_len=-10)
#     stim_ict = ti['stim_ict'][df_stimulus['stim_idx'].to_numpy() - 1]
#     df_stimulus['stim_ict'] = stim_ict
#
#     peak_mean, baseline_mean, baseline_std = respvec.compute_interval_response_vectors(C_dec,
#                                                                                        timestamps,
#                                                                                        df_stimulus['stim_ict'],
#                                                                                        **respvec_ops)
#     peak_amp = peak_mean - baseline_mean
#
#     # %% PCA, over odors
#     plt.close('all')
#     run_odor_pca = True
#     if run_odor_pca:
#         pca = decomposition.PCA(n_components=10)
#         Y = pca.fit_transform(preprocessing.normalize(X_peaks, norm='l2'))
#
#         with sns.axes_style("darkgrid"):
#             fig, axarr = plt.subplots(1, 2, tight_layout=True, figsize=(8, 4))
#             sns.scatterplot(x=Y[:, 0], y=Y[:, 1], hue=nt_meta_data.df_stimulus.stimname.to_list(),
#                             palette=_main.clr_dict,
#                             ax=axarr[1], alpha=1)
#             sns.lineplot(x=np.arange(pca.n_components), y=pca.explained_variance_ratio_, marker='o', ax=axarr[0])
#             fig.suptitle(title_str + """\n odor pca on peak_amp(C_dec)""", fontsize=10)
#             plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
#             plt.show()
#
#     # %% Run kenyon cell PCA
#     run_kc_pca = False
#     if run_kc_pca:
#         pca_kcs = decomposition.PCA()
#         Y_kcs = pca_kcs.fit_transform(peak_amp[good_cells_bool & bin_idx, :])
#         normed_comp = preprocessing.normalize(pca_kcs.components_, 'l2')
#
#         # plot first 2 components
#         with sns.axes_style("darkgrid"):
#             fig, axarr = plt.subplots(1, 2, tight_layout=True)
#             sns.scatterplot(x=Y_kcs[:, 0], y=Y_kcs[:, 1], ax=axarr[1])
#             sns.lineplot(x=np.arange(pca_kcs.n_components_) + 1, y=pca_kcs.explained_variance_ratio_, marker='o',
#                          ax=axarr[0])
#             plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
#             plt.show()
#
#         # plot KC components
#         df_plot = pd.DataFrame(normed_comp.T)
#         df_plot['stimname'] = df_stimulus.stimname
#         df_plot['stim_idx'] = df_stimulus.stim_idx
#         df_plot.sort_values(by=['stimname', 'stim_idx'], inplace=True)
#         df_plot.set_index(['stimname', 'stim_idx'], inplace=True)
#
#         fig, ax = plt.subplots(figsize=(7, 5))
#         ax = sns.heatmap(df_plot.transpose(), vmin=-.8, vmax=.8, cmap='RdBu_r')
#         ax.set_title('KC pca components, L2-normed')
#         plt.subplots_adjust(bottom=0.3, top=0.9, fontsize=10)
#         plt.show()
#
#     # %% plot ecdf of peak amplitudes for each stimulus
#     plt.close('all')
#     plot_peak_ecdf = False
#     if plot_peak_ecdf:
#         df_plot = pd.DataFrame(peak_mean[good_cells, :].transpose())
#         df_plot['stimname'] = df_stimulus['stimname']
#         df_plot['stim_idx'] = df_stimulus['stim_idx']
#         df_plot['reps'] = df_stimulus['reps']
#         df_plot = pd.melt(df_plot, id_vars=['stimname', 'stim_idx', 'reps'])
#
#         # ax = sns.ecdfplot(data=df_plot, x='value', hue='stimname')
#         with sns.axes_style('darkgrid'):
#             g = sns.displot(
#                 data=df_plot, x='value', col="stimname", hue='reps',
#                 palette=sns.dark_palette("#69d", reverse=False, as_cmap=True),
#                 kind="ecdf", height=4, aspect=.7)
#             g.fig.suptitle(title_str)
#             plt.subplots_adjust(bottom=0.15, top=0.8, hspace=0.1)
#             for ax in g.axes.flat:
#                 ax.set_ylim([0.7, 1.1])
#             plt.show()
#
#         # fig, ax = plt.subplots(1, 1)
#         # for i in range(df_stimulus.shape[0]):
#         #     ax.hist(peak_amp[good_cells, i], np.arange(0, 25, .1), density=True, histtype='step', cumulative=True,
#         #             color=clr_dict[df_stimulus.loc[i, 'stimname']], label=df_stimulus.loc[i, 'stimname'])
#         #
#         # plt.show()
#
#     # %% SVM: A vs B trained, use all 1-6ol and EP trials
#
#     plt.close('all')
#     include_odors = ['1-6ol', 'EP']
#     mask = df_stimulus.stimname.isin(include_odors).to_numpy()
#     mask_labels = df_stimulus.loc[mask, 'stimname'].to_list()
#
#     le = preprocessing.LabelEncoder()
#     le.fit_transform(mask_labels)
#
#     X = peak_amp[good_cells, :].T
#     X = X[mask, :]
#     X = preprocessing.normalize(X, )
#     y = le.transform(mask_labels)
#
#     clf = svm.SVC(kernel='linear')
#     clf.fit(X, y)
#     dec = clf.decision_function(preprocessing.normalize(peak_amp[good_cells, :].T))
#
#     df_svm = copy.deepcopy(df_stimulus)
#     df_svm['dec'] = dec
#     tdec = clf.decision_function(preprocessing.normalize(C_dec[good_cells, :].T))
#
#     save_svm_output = True
#     if save_svm_output:
#         np.save(f.with_name('svm_trial_results.npy'),
#                 dict(X_train=X, y_train=y, le=le, df_svm=df_svm, tdec=tdec),
#                 allow_pickle=True)
#
#     # save peak_amp svm trial - barplots
#     save_svm_barplot = True
#     if save_svm_barplot:
#         g = sns.catplot(data=df_svm, x='reps', y='dec', hue='stimname', col='stimname', kind='bar', dodge=False,
#                         height=5, aspect=0.2, legend_out=True)
#         g.set_titles("{col_name}")
#         g.set_axis_labels("rep", "decision fcn")
#         g.fig.suptitle(
#             title_str + """\nDecision function for peak amp. of deconvolved traces\n1-hexanol vs. ethyl propionate""",
#             fontsize=10)
#         plt.subplots_adjust(wspace=0.2, top=0.8)
#         sns.despine(trim=True)
#         plt.show()
#         g.fig.savefig(f.with_name('svm_trial.png'))
#
#     # save full timeseries, projected onto svm hyperplane
#     save_svm_timeseries = True
#     if save_svm_timeseries:
#         fig, ax = plt.subplots(figsize=(6, 3))
#         ax.fill_between(np.arange(tdec.size), tdec, np.zeros(tdec.shape), color='k')
#         fig.suptitle(
#             title_str + """\nDecision function for full deconvolved traces\nSVM trained on 1-hexanol vs. ethyl
#             propionate""",
#             fontsize=10)
#         plt.subplots_adjust(wspace=0.2, top=0.8)
#         plt.show()
#         fig.savefig(f.with_name('svm_timeseries.png'))
#
#     return ainfo, nt_meta_data, df_svm

# %%
#
#     skf = model_selection.StratifiedKFold(n_splits=3, shuffle=True)
#
#     for train_index, test_index in skf.split(X, y):
#         X_train, X_test = X[train_index], X[test_index]
#         y_train, y_test = y[train_index], y[test_index]
#         clf = svm.SVC(kernel='linear')
#         clf.fit(X, y)
#         dec = clf.decision_function(peak_amp[good_cells, :].T)
#
# #%%
#     X_train, X_test, y_train, y_test = train_test_split(X, y,
#                                                         test_size=0.5,
#                                                         random_state=42,
#                                                         shuffle=True)
#
#     # 1-6ol maps to 1
#
#     clf_hex = svm.SVC(kernel='linear')
#     clf_hex.fit(X_train, np.isin(y_train, [1, 2]))
#     decision_hex = clf_hex.decision_function(zdfof[iscell, :].transpose())
#
#     # %%
#     fig, axarr = plt.subplots(len(le.classes_), 1, figsize=(12, 8))
#
#     for ii, ax in enumerate(axarr.flat):
#         ax.plot(timestamps, decision_hex)
#         ax.set_title(le.classes_[ii])
#
#         clr_dict = {'1-5ol': 'tab:blue',
#                     '1-6ol': 'tab:green',
#                     'EP': 'tab:red',
#                     '1-6ol+EP': 'tab:purple',
#                     'pfo': 'gray'}
#
#         for row in df_stimulus.itertuples():
#             ax.axvspan(row.stim_ict + peak_start, row.stim_ict + peak_len, facecolor=clr_dict[row.stimname],
#             alpha=0.5)
#             if ii == 0:
#                 ax.text(row.stim_ict, 1.2, row.stimname, ha="left", size=8, rotation=45)
#         ax.axhline(0)
#
#     plt.show()
#     # %%
#
#     # EP maps to 3
#     clf_ep = svm.SVC(kernel='linear')
#     clf_ep.fit(X_train, np.isin(y_train, [2, 3]))
#     decision_ep = clf_ep.decision_function(zdfof[iscell, :].transpose())
#
#     # %%
#     fig, axarr = plt.subplots(len(le.classes_), 1, figsize=(12, 12))
#
#     for ii, ax in enumerate(axarr.flat):
#         ax.plot(timestamps, decision_ep)
#         ax.set_title(le.classes_[ii])
#
#         clr_dict = {'1-5ol': 'tab:blue',
#                     '1-6ol': 'tab:green',
#                     'EP': 'tab:red',
#                     '1-6ol+EP': 'tab:purple',
#                     'pfo': 'gray'}
#
#         for row in df_stimulus.itertuples():
#             ax.axvspan(row.stim_ict + peak_start, row.stim_ict + peak_len, facecolor=clr_dict[row.stimname],
#             alpha=0.5)
#             if ii == 0:
#                 ax.text(row.stim_ict, 1.2, row.stimname, ha="left", size=8, rotation=45)
#         ax.axhline(0)
#
#     plt.show()
#     # %%
#
#     clf = svm.SVC(kernel='linear', decision_function_shape='ovr', probability=True)
#     clf.fit(X_train, y_train == 1)
#
#     # compute decision function
#     Z = clf.decision_function(zdfof[iscell, :].transpose())
#     Zprob = clf.predict_proba(zdfof[iscell, :].transpose())
#
#     score = clf.score(X_test, y_test)
#     y_pred = clf.predict(X_test)
#     cm = confusion_matrix(y_test, y_pred, labels=clf.classes_, normalize='true')
#
#     # %%
#     Z = [None] * T
#     for ii in range(T):
#         Z[ii] = clf.decision_function(zdfof[iscell, ii].reshape(1, -1))
#
#     # %%
#
#     # %% load tensor information
#     tensor_raw = np.load(f.with_name('TENSOR_RAW.npz'))
#     tensor_F = tensor_raw['tensor_F']
#     tensor_Fneu = tensor_raw['tensor_Fneu']
#     trial_ts = tensor_raw['trial_ts']
#
#     tensor_Fc = tensor_F - tensor_Fneu
#     # %%
#     ustim = df_stimulus.stimname.unique().categories
#
#     # calculate means across odors
#     mean_odors = []
#     for iodor, odor in enumerate(ustim):
#         mean_odors.append(np.mean(tensor_Fc[:, df_stimulus.stimname == odor, :], axis=1))
#
#     mean_odor_tensor = np.stack(mean_odors)
#     mean_odor = np.hstack(mean_odors)
#
#     x = np.arange(mean_odor.shape[0])
#     ys = mean_odors + x[:mean_odor.shape[0], np.newaxis]
# %%

# %%#
# -----------------------------------------------------------------------------------------------------------------------
#  Set ordering of odors
#  Use to create a categorical dtype (ordered) in pandas. This can be used to set the dtype in df_stimulus['stimname'],
#  for easy ordering.

#  To order by stimname, do df_stimulus.sort_values(by='stimname')
#
# ocat = pd.CategoricalDtype(categories=['pfo', '1-5ol', 'VA', '1-6ol', 'EP', '1-6ol+EP'], ordered=True)
#
# color_dict = {'1-6ol': 'tab:blue',
#               '1-5ol': 'tab:green',
#               'EP': 'tab:red',
#               '1-6ol+EP': 'tab:purple',
#               'VA': 'tab:brown',
#               'pfo': 'tab:gray'}
#
# # namedtuple, structure for holding metadata
# MetaData = namedtuple("MetaData",
#                       ["df_stimulus",
#                        "df_frametimeinfo",
#                        "timestamps",
#                        "session_info",
#                        "xml_meta",
#                        "sync_meta",
#                        "meta_yaml",
#                        "fs", "blk_conc", "rear_type"
#                        ])
#
# Hallem = hallem.HallemSet()
#
#
# def get_title_str(ainfo_, nt_meta_data_):
#     """
#     Constructs title for plotting/analysis
#
#     :param ainfo_: info returned by dirutils.parse_dataset_info('.../suite2p/**/stat.npy')
#     :type ainfo_: namedtuple AnalysisSet in ryim.dirutils
#     :param nt_meta_data_: holds timing, stimulus, recording info
#     :type nt_meta_data_: namedtuple MetaData
#     :return: f-string, with relevant movie/block information
#     :rtype: str
#
#     """
#     title_str = f"""{ainfo_.flydate} fly {ainfo_.flynum}, {ainfo_.mov}"""
#
#     if ainfo_.is_block:
#         title_str = f"""{title_str} blk={ainfo_.blk}\n"""
#
#     title_str = title_str + f"""reared: {nt_meta_data_.rear_type}    |   conc={nt_meta_data_.blk_conc}"""
#     return title_str
#
#
# def postprocess_metadata(_nt_meta_data, ocat=ocat):
#     """
#     fix some of the metadata stuff
#     :param _nt_meta_data:
#     :param ocat:
#     :return:
#     """
#
#     _ti = _nt_meta_data.sync_meta
#     _df_stimulus = _nt_meta_data.df_stimulus
#     _df_frametimeinfo = _nt_meta_data.df_frametimeinfo
#
#     _df_stimulus['blk_conc'] = min(_df_stimulus.loc[:, ['conc_a', 'conc_b']].mode(axis=0).min(axis=0).tolist())
#     _df_stimulus['stim_ict'] = _ti['stim_ict'][_df_stimulus['stim_idx'].to_numpy() - 1]
#     _df_stimulus['scope_pulse_idx'] = np.interp(_df_stimulus['stim_ict'],
#                                                 _df_frametimeinfo['frame_times'],
#                                                 _df_frametimeinfo['scope_pulse_idx']).astype('int')
#     _df_stimulus['stimname'] = _df_stimulus['stimname'].astype(ocat)
#
#     _df_stimulus['reps'] = ''
#     for ustim in _df_stimulus['stimname'].unique():
#         mask = _df_stimulus['stimname'] == ustim
#         _df_stimulus.loc[mask, 'reps'] = np.arange(np.count_nonzero(mask*1))
#
#     _session_info = {key: value for key, value in _nt_meta_data.session_info.to_dict().items() if key.isidentifier()}
#
#     _fs = trace.compute_fs(_nt_meta_data.timestamps)
#     _blk_conc = _df_stimulus.blk_conc.unique()[0]
#     _rear_type = _session_info['rearing']
#
#     return MetaData(_df_stimulus,
#                     _df_frametimeinfo,
#                     _nt_meta_data.timestamps,
#                     _session_info,
#                     _nt_meta_data.xml_meta,
#                     _nt_meta_data.sync_meta,
#                     _nt_meta_data.meta_yaml,
#                     _fs,
#                     _blk_conc,
#                     _rear_type
#                     )
#
#
# def make_label_for_svm(ts, stim_ict, stim_idx, peak_start, peak_len):
#     """
#     Fill a vector (like timestamps) w/ stim_idx at specified intervals w/ after stim_ict.
#
#     :param ts: timestamps
#     :type ts: np.ndarray, 1D
#     :param stim_ict: time that stimulus turns on; df_stimulus.stim_ict
#     :type stim_ict: np.ndarray
#     :param stim_idx: 1-indexed stimulus # index (ex: df_stimulus.stim_idx)
#     :type stim_idx:  np.ndarray
#     :param peak_start: time (in seconds) after stim_ict the peak window starts
#     :type peak_start: float
#     :param peak_len: time (in seconds) after stim_ict + peak_start, how long the peak window lasts for
#     :type peak_len: float
#     :return: labeltimestamp
#     :rtype:  np.array
#
#     """
#     label_ts = np.zeros_like(ts)
#
#     for s_ict, s_idx in zip(stim_ict, stim_idx):
#         t0 = s_ict + peak_start
#         t1 = s_ict + peak_start + peak_len
#         idx = (ts >= t0) & (ts < t1)
#         label_ts[idx] = s_idx
#     return label_ts
#
#
# def svm_a_vs_b(ts, stim_a_ict, stim_b_ict, peak_start, peak_len):
#     """ One classifier, trained to discriminate between odor a and odor b."""
#     pass
#
#
#
# def draw_oasis_results(y_, c_, sp_, x_=None, sn_=None, ax_=None):
#     """
#     Plots fluorescence activity and results of calcium deconvolution.
#
#     :param y_: fluorescence trace, raw
#     :param c_: deconvolved trace
#     :param sp_: spikes
#     :param x_: timestamps for y_, c_, sp_
#     :param sn_: SNR, scalar
#     :param ax_: axes
#     :return: axes w/ plots
#     """
#
#     if x_ is None:
#         x_ = np.arange(y.size)
#
#     if ax_ is None:
#         ax_ = sns.lineplot(x=x_, y=y_, color='k', alpha=0.3)
#     else:
#         sns.lineplot(x=x_, y=y_, color='k', alpha=0.3, ax=ax_)
#
#     sns.lineplot(x=x_, y=c_,  color='k', alpha=1, ax=ax_)
#     sns.rugplot(x=x_[sp_ > 0], height=0.1, color='r', ax=ax_)
#
#     if sn_ is not None:
#         ax_.axhline(sn_, linestyle='--', alpha=0.5)
#
#     return ax_
#
#
# def load_metadata(f, method='statfile'):
#     """Helper func"""
#     if method == 'statfile':
#         ainfo = dirutils.parse_dataset_info(f)
#
#     elif method == 'AnalysisSet':
#         ainfo = f
#
#     print(ainfo)
#     # loading metadata
#     if ainfo.is_block:
#         meta_data_loader = analysis.BlockMetaDataLoader.init_from_statfile(f)
#     else:
#         meta_data_loader = analysis.MovieMetaDataLoader.init_from_statfile(f)
#     return meta_data_loader.load()
#
#
# def load_suite2p_data(f, method='statfile'):
#     """Helper func for eaasy loading of data"""
#     if method == 'statfile':
#         ainfo = dirutils.parse_dataset_info(f)
#     elif method == 'AnalysisSet':
#         ainfo = f
#
#     # loading suite2p output data
#     if f.parent.name == 'combined':
#         suite2p_loader = analysis.Suite2pCombinedLoader(statfile=f)
#     elif 'plane' in f.parent.name:
#         suite2p_loader = analysis.Suite2pPlane0Loader(statfile=f)
#
#     nt_s2p_data_ = suite2p_loader.load()
#     return nt_s2p_data_
#
#
# def load_all(f):
#     ainfo = dirutils.parse_dataset_info(f)
#     nt_meta_data = postprocess_metadata(load_metadata(f))
#     nt_s2p_data = load_suite2p_data(f)
#     title_str = get_title_str(ainfo, nt_meta_data)
#     return nt_meta_data, nt_s2p_data, title_str
