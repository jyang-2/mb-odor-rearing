#  STEP 1. Setup

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import nitime
# import nitime.timeseries as ts
# import nitime.analysis as nta
# import nitime.viz as viz
import importlib as imp
import json
from scipy.stats import zscore
import json
import pprint as pp

import rasterplot
import copy
import plotly.express as px

import ryim.chonk as chonk
from ryim.dirutils import *
from rearing.ca import respvec, trace, vis

from suite2p.extraction import dcnv
import numpy as np
from rearing.ca import respvec, trace, vis
from scipy import signal

import plotly.express as px
import plotly.graph_objects as go
import plotly.subplots as sp
import plotly.io as pio

pio.renderers.default = "browser"
from rearing import tidy

# %%
import matplotlib as mpl

mpl.rcParams.update({
    'axes.spines.left': True,
    'axes.spines.bottom': True,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'legend.frameon': False,
    'ytick.major.left': True,
    'figure.subplot.wspace': .02,
    'figure.subplot.hspace': .02,
})


# %%

def print_dir_info(file):
    print(file)
    print(f"\tflydate: {parse_parent_date(file).name}")
    print(f"\tflynum: {parse_parent_fly(file).name}")
    print(f"\tmovie: {parse_parent_movie(file).name}")
    print(f"\tis in block folder: {in_block_dir(file)}")
    print(f"\tblock: {parse_parent_block(file).name}")
    return


# %% Select list of datasets (look for embedding.npy)
"""
Choosing one of these files also selects a date, fly, and movie to run.
"""

DATA_DIR = Path('./data/processed_data')
# files = sorted(list(DATA_DIR.rglob(f"2021-08-30/**/combined/embedding.npy")))
files = sorted(list(DATA_DIR.rglob(f"2021-08-30/**/embedding.npy")))

list_files = True
if list_files:
    for i, f in enumerate(files):
        print(f"{i} ----- {f}")

# %% Select one file to run

# f = Path("../data/processed_data/2021-08-30/2/movie_001/1/suite2p/combined/embedding.npy")
f = files[0]

# Print dataset information
is_blk = in_block_dir(f)
print_dir_info(f)

# %%  Load data for selected file.
"""

** Note: rear_type calculated here **

Parse the filepath to ** embedding.npy ** to extract info about  the  dataset, including rearing condition.
	- `flydate`
	- `flynum`
	- `mov`
	- `blk`

Load analysis outputs  from the suite2p folder.
	- `Fcell` -from F.npy
	- `Fneu` - from Fneu.npy
	- `iscell` - from iscell.npy (first column = 0 or 1, second column = probability of being a cell)
	- `stat` - from stat.npy
	- `embedding` - from embedding.npy

Depending on if the dataset 's movie was split into chunks (sub-movies) for analysis, 
load `df_frametimeinfo`, `df_stimulus`, and `timestamps` from the chunk directory (`blk_dir` )
or the movie directory (`blk_dir`)
"""

if is_blk:
    blk_dir = parse_parent_block(f)
    mov_dir = parse_parent_movie(f)

    df_stimulus = pd.read_csv(blk_dir.joinpath('metadata.csv'))
    df_frametimeinfo = pd.read_csv(blk_dir.joinpath('frametimeinfo.csv'))
    timestamps = np.load(blk_dir.joinpath('timestamps.npy'))

    session_info = pd.read_excel(mov_dir.joinpath('metadata.xlsx'), sheet_name='sessioninfo').squeeze()
    rear_type = f"{session_info['rearing']} reared"

    flydate = parse_parent_date(f).name
    flynum = parse_parent_fly(f).name
    mov = parse_parent_movie(f).name
    blk = parse_parent_block(f).name
    title_str = f"{flydate} fly {flynum}, {mov} block {blk}, {rear_type}"

    # load files from suite2p folder
    s2p_dir = blk_dir.joinpath('suite2p')
    Fcell = np.load(f.with_name('F.npy'))
    Fneu = np.load(f.with_name('Fneu.npy'))
    iscell = np.load(f.with_name('iscell.npy'))
    stat = np.load(f.with_name('stat.npy'), allow_pickle=True)
    embedding = np.load(f.with_name('embedding.npy'), allow_pickle=True).item()
else:  # is_blk: True
    mov_dir = parse_parent_movie(f)

    df_stimulus = pd.read_excel(mov_dir.joinpath('metadata.xlsx'), sheet_name='stimulus')
    df_frametimeinfo = pd.read_csv(mov_dir.joinpath('frametimeinfo.csv'))
    timestamps = np.load(mov_dir.joinpath('timestamps.npy'))

    session_info = pd.read_excel(mov_dir.joinpath('metadata.xlsx'), sheet_name='sessioninfo').squeeze()
    rear_type = f"{session_info['rearing']} reared"

    flydate = parse_parent_date(f).name
    flynum = parse_parent_fly(f).name
    mov = parse_parent_movie(f).name
    title_str = f"{flydate} fly {flynum}, {mov}, {rear_type}"

    # load files from suite2p folder
    s2p_dir = mov_dir.joinpath('suite2p')
    Fcell = np.load(f.with_name('F.npy'))
    Fneu = np.load(f.with_name('Fneu.npy'))
    iscell = np.load(f.with_name('iscell.npy'))
    stat = np.load(f.with_name('stat.npy'), allow_pickle=True)
    embedding = np.load(f.with_name('embedding.npy'), allow_pickle=True).item()

# calculate fluorescence
F = Fcell - 0.7 * Fneu
F = F[iscell[:, 0].astype(np.bool_), :]
n_cells = F.shape[0]

# normalize/preprocess neuropil-corrected fluorescence traces
sp = copy.deepcopy(F)
spn = rasterplot.norm_traces(sp)
spp = rasterplot.prep_to_plot(spn)

# compute block/chunk concentration, from df_stimulus
blk_conc = min(df_stimulus.loc[:, ['conc_a', 'conc_b']].mode(axis=0).min(axis=0).tolist())

print(f"\n{title_str}")
print(f"block conc = {blk_conc}")
# print(f"\ndf_stimulus:")
print(f"\n_______df_stimulus__________:")
print(df_stimulus)

print(f"_______df_frametimeinfo__________:")
print(df_frametimeinfo.head())

print(f"\nsize of timestamps array: {timestamps.shape}")
print(f"\nsize of df_frametimeinfo dataframe: {df_frametimeinfo.shape}")

# %% Correct df_frametimeinfo
"""
If df_frametimeinfo contains timestamps for multiple planes, drop everything but plane=0 timestamps,
and check that the dimensions match timestamps
"""

if 'plane' in df_frametimeinfo.columns:
    df_frametimeinfo = df_frametimeinfo.loc[df_frametimeinfo['plane'] == 0]
    df_frametimeinfo.reset_index(drop=True, inplace=True)  #

print(f"\nsize of timestamps array: {timestamps.shape}")
print(f"\nsize of df_frametimeinfo dataframe: {df_frametimeinfo.shape}")

print(df_frametimeinfo.head())

# %%
"""
From sync_meta.json file, load into dict ti, w/ fields

- block_fct : start time for each acquisition
- block_ict: end time for acquisition
- n_blocks: # of acquisitions
- n_stim: total # of stimuli presented
- stim_ict: stimulus onset times
- stim_fct: stimulus offset times

Note: this is always loaded from the movie folder, not the block folder. 

If you are analyzing a block, then go to the parent folder to load, and note
 that it will include stimuli info for different blocks as well.
"""
# movie folder should have a sync_meta.json file, containing stimulus and acquisition onset/offset times
# load into dict ti (trial info)
with open(mov_dir.joinpath('sync_meta.json')) as file:
    ti = json.load(file)

# make everything into a 1d numpy array
for k, v in ti.items():
    if isinstance(v, list):
        ti[k] = np.array(v)

pp.pprint(ti)

# %%

# calculate timing of stimuli, using ti['stim_ict']
stim_ict = ti['stim_ict'][df_stimulus['stim_idx'].to_numpy() - 1]
df_stimulus['stim_ict'] = stim_ict

# show df_stimulus, with new column 'stim_ict'
print(df_stimulus)

# %% Set ordering of odors
"""
Set ordering of odors
Use to create a categorical dtype (ordered) in pandas. This can be used to set the dtype in df_stimulus['stimname'], for easy ordering.

To order by stimname, do df_stimulus.sort_values(by='stimname')
"""

# make categorical datatype
ocat = pd.CategoricalDtype(categories=['pfo', '1-5ol', 'VA', '1-6ol', 'EP', '1-6ol+EP'], ordered=True)
df_stimulus['stimname'] = df_stimulus['stimname'].astype(ocat)

print("df_stimulus, sorted by column 'stimname'")
print(df_stimulus.sort_values(by='stimname'))

# %%
s2p_preproc_ops = dict()
tau = 20.0  # timescale of indicator
fs = trace.compute_fs(timestamps)  # sampling rate in Hz
neucoeff = 0.7  # neuropil coefficient
# for computing and subtracting baseline
baseline = 'maximin'  # take the running max of the running min after smoothing with gaussian
sig_baseline = 10.0  # in bins, standard deviation of gaussian with which to smooth
win_baseline = 60.0  # in seconds, window in which to compute max/min filters

# split traces
blocked_traces = trace.split_traces_to_blocks(F, df_frametimeinfo.scope_pulse_idx)

# lowpass filter
filt_ops = dict(fs=fs, cutoff=0.5, trans_width=0.05, numtaps=230)
b = signal.remez(filt_ops['numtaps'], [0, filt_ops['cutoff'], filt_ops['cutoff'] + filt_ops['trans_width'],
                                       0.5 * filt_ops['fs']], [1, 0], Hz=filt_ops['fs'])
filt_blocked_traces = [signal.filtfilt(b, 1, y) for y in blocked_traces]

dfof_method = 'trace'

# suite2p preprocessing
if dfof_method == 'suite2p':
    blocked_dfof, blocked_df, blocked_baseline = zip(*[trace.suite2p_extract_dfof(
        F=item,
        baseline='constant_prctile',
        win_baseline=90,
        sig_baseline=5,
        fs=fs,
        prctile_baseline=10) for item in filt_blocked_traces])

# rolling percentile window
if dfof_method == 'trace':
    blocked_dfof, blocked_df, blocked_baseline = zip(*[trace.extract_dfof(item,
                                                                          fs=fs,
                                                                          win=60,
                                                                          percentile=10)
                                                       for item in blocked_traces])

blocked_zdfof = [zscore(item, axis=1) for item in filt_blocked_traces]

# combine lists of traces
dfof = np.hstack(tuple(blocked_dfof)).astype('float32')
zdfof = np.hstack(tuple(blocked_zdfof)).astype('float32')

# %% compute peak response amplitudes

peak_start = 1
peak_len = 5

baseline_len = -10

n_trials = df_stimulus.shape[0]

peak_mean = np.empty((n_cells, n_trials))
baseline_mean = np.empty((n_cells, n_trials))
baseline_std = np.empty((n_cells, n_trials))

for i, stim_on in enumerate(df_stimulus['stim_ict']):
    peak_idx = (timestamps >= (stim_on + peak_start)) \
               & (timestamps < (stim_on + peak_start + peak_len))
    peak_mean[:, i] = np.mean(zdfof[:, peak_idx], axis=1)

    baseline_idx = (timestamps >= (stim_on + baseline_len)) & (timestamps < stim_on)
    baseline_mean[:, i] = np.mean(zdfof[:, baseline_idx], axis=1)
    baseline_std[:, i] = np.std(zdfof[:, baseline_idx], axis=1)

peak_amp = peak_mean - baseline_mean
peak_nstd = peak_amp / baseline_std

# compute responders
bin_response = respvec.get_binresponse(peak_amp, baseline_std, n_std=3, thr=0.5)
bin_idx = np.any(bin_response, axis=1)

dist_all, df_peakresp = respvec.compute_corrmats(peak_amp, df_stimulus)
dist_resp, df_peakresp = respvec.compute_corrmats(peak_amp[bin_idx, :], df_stimulus)

# %% plot correlation matrices

# fig corrmat_all: plot correlation matrices, with all cells included
# ----------------------------------------------------------------------

fig_corrmat_all, axarr, haxarr = vis.corrmats(dist_all)
fig_corrmat_all.suptitle(f"{title_str}\nconc={blk_conc}, reared={rear_type}\nall cells")

# fig_corrmat_resp: plot correlation matrices, with only responders included
# ----------------------------------------------------------------------

fig_corrmat_resp, axarr, haxarr = vis.corrmats(dist_resp)
fig_corrmat_resp.suptitle(f"{title_str}\nconc={blk_conc}, reared={rear_type}\n responders only")
plt.show()

# save corrmat plots
fig_corrmat_resp.savefig(f.with_name('corrmat_resp.pdf'))
fig_corrmat_all.savefig(f.with_name('corrmat_all.pdf'))
# %%

tidy_pearson_all = tidy.tidy_triu_df_corrmat(dist_all['pearson'])
print(tidy_pearson_all.head())

# %% plot traces


acq_idx = np.interp(df_stimulus['stim_ict'],
                    df_frametimeinfo['frame_times'],
                    df_frametimeinfo['scope_pulse_idx']).astype('int')
df_stimulus['scope_pulse_idx'] = acq_idx

cids = list(range(20))

# plot traces in subplot array
fig, axarr = vis.blocked_stackplot(blocked_zdfof, timestamps, cids=cids,
                                   df_stimulus=df_stimulus, df_frametimeinfo=df_frametimeinfo,
                                   peak_start=peak_start, peak_len=peak_len)

fig.suptitle(f"{title_str}\nconc={blk_conc}, reared={rear_type}")
plt.show()
# %%
skip = True
if ~skip:
    fig, axarr = plt.subplots(len(cids), nacq, figsize=(12, 9), sharey=True, sharex='col')

    scope_pulse_idx = df_frametimeinfo['scope_pulse_idx'] - 1
    acq_idx = np.interp(df_stimulus['stim_ict'],
                        df_frametimeinfo['frame_times'],
                        df_frametimeinfo['scope_pulse_idx']).astype('int')
    for cid in cids:
        for i in range(nacq):
            ts = timestamps[scope_pulse_idx == i]
            print(ts[0])
            print(ts[-1])
            axarr[cid, i].plot(ts, blocked_traces[i][cid, :], label='blocked_traces')
            axarr[cid, i].plot(ts, filt_blocked_traces[i][cid, :], label='filt_blocked_traces')
            axarr[cid, i].plot(ts, blocked_baseline[i][cid, :], label='blocked_baseline')

    for ax in axarr.flat:
        ax.label_outer()
    plt.ylabel('dfof', fontsize=9)

    for idx, col, ostim in zip(df_stimulus['stim_ict'], acq_idx - 1, df_stimulus['stimname']):  #
        print(f"{idx},{col},{ostim}")
        axarr[0, col].axvline(idx, label=ostim)
    plt.legend(bbox_to_anchor=(1, 1),
               bbox_transform=plt.gcf().transFigure)

# %% plotly trace plot

fig = vis.dfof_facet_plot(dfof, df_frametimeinfo, cids=list(range(10)), yrange=(-1.0, 8.0))

# print(f"{idx}, {col}, {ostim}")
# fig.add_vline(x=idx,
#               row='all', col=col,
#               line_dash="dot",
#               annotation_text=ostim)
fig.show()

# %%
import seaborn as sns

sns.relplot(
    data=dfof[:10, :],
    x="time", y="firing_rate",
    hue="coherence", size="choice", col="align",
    kind="line", size_order=["T1", "T2"], palette=palette,
    height=5, aspect=.75, facet_kws=dict(sharex=False),
)
# %% md

# Plotting
## Draw rastermaps
"""
1 plot per odor delivered, with stimulus onsets indicated & labeled.

Plots should be identical to what you see when you run the suite2p GUI rastermap module.
"""


# %%

# suite2p rastermap plotting function
def suite2p_display(S, isort=None, nbin=50, figure_args={}, imshow_args={}):
    """
		Reproduces suite2p rastermap visualization.

		S = cell responses (# cells x timepoints)
	"""
    if isort is not None:
        S = S[isort, :]

    Sfilt = rasterplot.running_average(S, nbin)
    Sfilt = zscore(Sfilt, axis=1)

    # fig1 = plt.figure(figsize=(16, 8), **figure_args)
    fig = px.imshow(Sfilt, zmin=-0.5, zmax=5, aspect='auto', color_continuous_scale='gray_r', origin='lower',
                    **imshow_args)

    plt.xlabel('time points')
    plt.ylabel('sorted neurons')
    return fig


# %%

# %%

if spp.shape[0] < 500:
    nbin = 1
elif spp.shape[0] > 1000:
    nbin = 30

# %%

# sp = copy.deepcopy(F)
# spn = rasterplot.norm_traces(sp)
# spp = rasterplot.prep_to_plot(spn)

# calculate sorting order of cells, from rastermap analysis
model_dict = vars(embedding)  # turn rastermap model into a dictionary for easy handling
isort1 = np.argsort(model_dict['embedding'][:, 0])  # calculate sort order from embedding values

odor_list = df_stimulus['stimname'].unique()
# blk_conc = min(df_stimulus.loc[:, ['conc_a', 'conc_b']].mode(axis=0).min(axis=0).tolist())

draw_plot = True  # whether or not to actually make the plots
for odor in odor_list:

    df = df_stimulus.loc[df_stimulus['stimname'] == odor]
    stim_idx_ = df['stim_idx'].to_numpy() - 1
    stim_ict = np.array(ti['stim_ict'])
    stim_ts = stim_ict[stim_idx_]

    # print(f"{odor}: {stim_ts}")

    if draw_plot:
        # ------------PLOT RASTERMAP----------------
        fig, ax = rasterplot.suite2p_display(spp, isort=isort1, nbin=nbin)
        ax.set_title(f"{title_str}\n# cells={sp.shape[0]}, odor={odor},  conc={blk_conc}")
        stim_ici = np.floor(np.interp(stim_ts, df_frametimeinfo.frame_times, df_frametimeinfo.index))
        for t in stim_ici:
            ax.axvline(t, linewidth=1, linestyle='--')
            ax.text(t, sp.shape[0] * 0.9, f" {odor}", rotation=90, verticalalignment='center',
                    horizontalalignment='right')
        fig.show()

# %% md

## Draw reordered rastermaps, w/ stimuli grouped

# %%



# %%
