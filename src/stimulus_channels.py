import suite2p
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import seaborn as sns
from pathlib import Path
import pandas as pd
import pprint
from drosolf import orns, pns
from scipy import signal
import scipy.io as sio
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore
from scipy.ndimage import filters
from scipy.cluster.hierarchy import dendrogram, linkage, fclusterdata
import importlib
import textwrap
from rastermap.mapping import Rastermap
from itertools import product

from analysis_respvec import make_sparsity_long
from rearing import hallem, io, meta
from rearing.ca import vis, trace, respvec
from sklearn import preprocessing
from datetime import datetime
import heapq
import timeit

pp = pprint.PrettyPrinter(indent=2, compact=True)
Hallem = hallem.HallemSet()

params = {'legend.fontsize': 8,
          'axes.labelsize': 8,
          'axes.titlesize': 10,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'font.family': ['sans-serif'],
          'figure.dpi': 300.0
          }
plt.rcParams.update(params)

# %% Set dataset info and directory locations

# def main():

expt = dict()
expt['date'] = '2021-06-26'
expt['fly'] = 1
expt['microscope'] = 'galvo-res'

# make directories
PRJ_DIR = Path.cwd()
FLY_DIR = PRJ_DIR.joinpath('data', 'processed_data', expt['date'], str(expt['fly']))

# load metadata yaml
metadata = io.load_yaml(FLY_DIR.joinpath("metadata.yaml"))
acquisitions = metadata['Acquisitions']['fly1']

# %% set movie/acquisition number and load data

expt['movie'] = 'movie_002'
iacq = [item['thorimage'] for item in acquisitions].index(expt['movie'])

# make movie directories
PROC_DIR = FLY_DIR.joinpath(expt['movie'])
S2P_DIR = PROC_DIR.joinpath('suite2p', 'plane0')

vals = suite2p.gui.io.load_files(S2P_DIR.joinpath('stat.npy'))
keys = ('stat', 'ops', 'Fcell', 'Fneu', 'Spks', 'iscell', 'probcell', 'redcell', 'probredcell', 'hasred')
s2p_dat = dict(zip(keys, vals))
# %%
# load metadata pandas dfs
df_sessioninfo = pd.read_excel(PROC_DIR.joinpath('metadata.xlsx'), sheet_name='sessioninfo')
df_stimulus = pd.read_excel(PROC_DIR.joinpath('metadata.xlsx'), sheet_name='stimulus')
df_stimulus = df_stimulus.rename(columns={'Unnamed: 0': 'stim_idx'})
df_stimulus['stim_idx'] += 1
df_stimulus['stimname'] = pd.Categorical(df_stimulus['stimname'], ordered=True,
                                         categories=['empty', 'pfo', '1-5ol', '1-6ol', '1-5ol+1-6ol', 'EP'])

df_stimtype = pd.read_excel(PROC_DIR.joinpath('metadata.xlsx'), sheet_name='stimtype')

# load timing, trial info
df_frametimeinfo = pd.read_csv(PROC_DIR.joinpath('frametimeinfo.csv'))
ti0 = sio.loadmat(PROC_DIR.joinpath('ti.mat'))
ti = {key: value.squeeze().tolist() for key, value in ti0.items() if isinstance(value, np.ndarray)}
fti = df_frametimeinfo.to_records()
stim_frames = np.interp(ti['stim_ict'], fti.frame_times, fti.index)
block_frames = np.interp(ti['block_ict'], fti.frame_times, fti.index)

# load suite2p files
ops0 = np.load(S2P_DIR.joinpath('ops.npy'), allow_pickle=True).item()
iscell = np.load(S2P_DIR.joinpath('iscell.npy'), allow_pickle=True)
embedding = np.load(S2P_DIR.joinpath('embedding.npy'), allow_pickle=True).item()

F = np.load(S2P_DIR.joinpath('F.npy'), allow_pickle=True)
Fneu = np.load(S2P_DIR.joinpath('Fneu.npy'), allow_pickle=True)

F = F[iscell[:, 0] == 1, :]
Fneu = Fneu[iscell[:, 0] == 1, :]
Fc = F - 0.7 * Fneu

K, T = Fc.shape
n_trials = ti['num_stim']
n_blocks = ti['num_blocks']
movie_conc = df_stimulus.loc[:, ['conc_a', 'conc_b']].min().min()

# %% Create lowpass filter

fs = ops0['fs']  # Sample rate, Hz
cutoff = .5  # Desired cutoff frequency, Hz
trans_width = .1  # Width of transition from pass band to stop band, Hz
numtaps = 203  # Size of the FIR filter.
b = signal.remez(numtaps, [0, cutoff, cutoff + trans_width, 0.5 * fs], [1, 0], Hz=fs)

# plot frequency response of filter
w, h = signal.freqz(b, [1], 2000)
# trace.plot_filter_response(fs, w, h, "Low-pass Filter")

# lowpass filter Fc, block by block
trial_traces = np.array_split(Fc, ti['num_blocks'], axis=1)
filt_trial_traces = [signal.filtfilt(b, 1, y) for y in trial_traces]

# calculate df/f from unfiltered traces
filter_traces = True
if filter_traces:
    # calculate df/f from lowpass filtered traces
    baseline = [trace.percentile_baseline(y, fs=6.1, win=45, percentile=20) for y in filt_trial_traces]
    df = [y - b for y, b in zip(filt_trial_traces, baseline)]
    dfof = [y / b for y, b in zip(df, baseline)]
    dfof_z = [trace.preprocess_for_rastermap(item, 'suite2p') for item in dfof]
    filtdat = {'baseline': baseline, 'df': df, 'dfof': dfof, 'dfof_z': dfof_z}
else:
    baseline = [trace.percentile_baseline(y, fs=6.1, win=45, percentile=10) for y in trial_traces]
    df = [y - b for y, b in zip(trial_traces, baseline)]
    dfof = [y / b for y, b in zip(df, baseline)]
    dfof_z = [trace.preprocess_for_rastermap(item, 'suite2p') for item in dfof]
    unfiltdat = {'baseline': baseline, 'df': df, 'dfof': dfof, 'dfof_z': dfof_z}

# %%

# Rastermap works the same as TSNE or UMAP
# create an instance of the Rastermap class
# class Rastermap:
"""rastermap embedding algorithm
    Rastermap first takes the top iPC's of the data, and then embeds them into
    n_X clusters. It returns upsampled cluster identities (n_X*upsamp).
    Clusters are also computed across Y (n_Y) and smoothed, to help with fitting.
    If n_components=2, n_X x n_X components are used.

    Parameters
    -----------
    n_components: int, optional (default: 1)
        dimension of the embedding space
    n_X: int, optional (default: 40)
        number of clusters in X
    n_Y: int, optional (default: 0)
        number of clusters in Y: will be used to smooth data before sorting in X
        if set to zero, no smoothing will occur
    iPC: nparray, int, optional (default: 0-399)
        which PCs to use during optimization
    nPC: int, optional (default: 400)
        if used, will take 0-nPC PCs during optimization
    init : initialization of algorithm (default: 'pca')
        can use 'pca', 'random', or a matrix n_samples x n_components
"""
# rastermap options
ops = {'n_components': 1, 'n_X': 40, 'alpha': 1., 'K': 1.,
       'nPC': 200, 'constraints': 2, 'annealing': True, 'init': 'pca',
       'start_time': 0, 'end_time': -1}

model = Rastermap(n_components=ops['n_components'], n_X=ops['n_X'], nPC=ops['nPC'],
                  init=ops['init'], alpha=ops['alpha'], K=ops['K'],
                  constraints=ops['constraints'], annealing=ops['annealing'])
dfof = filtdat['dfof']
dfof = np.hstack(tuple(dfof)).astype('float32')
model.fit(dfof)
proc = {'model': model, 'ops': ops, 'activity_mode': 'filtered dfof', 'sp': dfof}
np.save(S2P_DIR.joinpath('rastermodel.npy'), proc)
# %% response vector calculation

dfof = filtdat['dfof']
#dfof = np.hstack(tuple(dfof)).astype('float32')

stim_list = np.arange(n_trials) + 1

peakresp_nmean, baseline_avg, baseline_std = respvec.compute_response_vectors(dfof, baseline_idx=fti.baseline_idx,
                                                                              peakresp_idx=fti.peakresp_idx,
                                                                              stim_list=stim_list,
                                                                              maxk=3,
                                                                              stack_numpy=True)
peakresp_amp = peakresp_nmean - baseline_avg
#%%


# %% fig 0 - sparsity plot over std_thr range

recalculate = True

if recalculate:
    n_std = np.arange(0, 20, 0.1)
    sparsity = respvec.compute_sparsity(n_std, peakresp_amp, baseline_std)

    # use mean empty sparsity to zero the matrix
    meansparsity_empty = np.mean(sparsity[:, np.where(stimname == "empty")[0]], axis=1, keepdims=True)
    sparsity0 = sparsity - meansparsity_empty

    # sparsity dataframe, long form
    df_sparsity = make_sparsity_long(sparsity,
                                     n_std=n_std,
                                     stim_idx=df_stimulus['stim_idx'].to_numpy(),
                                     stimname=df_stimulus['stimname'].to_numpy())
    df_sparsity['zeroed'] = False

    # sparsity dataframe w/ empty subtracted, long form
    df_sparsity0 = make_sparsity_long(sparsity0,
                                     n_std=n_std,
                                     stim_idx=df_stimulus['stim_idx'].to_numpy(),
                                     stimname=df_stimulus['stimname'].to_numpy())
    df_sparsity0['zeroed'] = True
# %%
sns.set_context('notebook')

df = df_sparsity.append(df_sparsity0)

g = sns.FacetGrid(df, col="stimname", row="zeroed", aspect=0.8, height=4)
g.map(sns.lineplot, "n_std", "sparsity")
g.map(plt.axhline, y=0.1, ls=":", c=".5")
g.map(plt.axhline, y=0.05, ls=":", c=".5")
g.fig.subplots_adjust(top=0.8)  # adjust the Figure in rp
g.fig.suptitle(f"date={expt['date']}, fly={expt['fly']}, movie={expt['movie']}, conc={movie_conc}\n"
               f"scope={expt['microscope']}, reared=1-5ol+1-6ol, DATA=peakresp(filtdat['dfof'])")
fig0 = g.fig
plt.show()

#%%
title_str = f"date={expt['date']}, fly={expt['fly']}, movie={expt['movie']}, conc={movie_conc}\n" \
            f"scope={expt['microscope']}, reared=1-5ol+1-6ol, DATA=peakresp(filtdat['dfof'])"

sns.set_style('darkgrid')
df = df_sparsity.append(df_sparsity0)
g = sns.FacetGrid(df, col="zeroed", hue='stimname', aspect=0.8, height=4, sharey=False)
g.map(sns.lineplot, "n_std", "sparsity")
#
# g = sns.relplot(data=df,
#                 x="n_std", y="sparsity",
#                 hue="stimname", kind='line')
g.map(plt.axhline, y=0.1, ls=":", c=".5")
g.map(plt.axhline, y=0.05, ls=":", c=".5")
g.fig.subplots_adjust(top=0.8)  # adjust the Figure in rp
g.fig.suptitle(title_str)
plt.show()

# %% fig 1 - correlations
# ----------------------------------------------------------------------
# calculate peakresp_amp distances (cosine, pearson, spearman, kendall)
# ----------------------------------------------------------------------
recalculate = True
if recalculate:
    df_stimulus['stimname'] = pd.Categorical(df_stimulus['stimname'], ordered=True,
                                             categories=['empty', 'pfo', '1-5ol', '1-6ol', '1-5ol+1-6ol', 'EP'])
    sort_idx = df_stimulus['stimname'].sort_values().index.to_list()

    df_peakresp = pd.DataFrame(peakresp_amp)
    df_peakresp.columns = pd.MultiIndex.from_frame(df_stimulus.loc[:, ['stim_idx', 'stimname']])
    df_peakresp = df_peakresp.iloc[:, sort_idx]

    dist = dict()
    dist['cosine'] = 1 - squareform(pdist(df_peakresp.transpose(), metric='cosine'))
    dist['pearson'] = df_peakresp.corr(method='pearson')
    dist['spearman'] = df_peakresp.corr(method='spearman')
    dist['kendall'] = df_peakresp.corr(method='kendall')

# ----------------------------
# seaborn plotting (heatmaps)
# -----------------------------
sns.set_context('paper')
ticklabels = [f"{item[1]} ({item[0]})" for item in df_peakresp.columns.to_list()]
odors = [item[1] for item in df_peakresp.columns.to_list()]
title_str = f"date={expt['date']}, fly={expt['fly']}, movie={expt['movie']}, conc={movie_conc}\n" \
            f"scope={expt['microscope']}, reared=1-5ol+1-6ol, DATA=peakresp(filtdat['dfof'])"
clim = [-1.0, 1.0]

fig1, axarr = plt.subplots(2, 2, figsize=(8, 8), constrained_layout=True)
axarr = axarr.flatten()
haxarr = [None] * 4

for i, m in enumerate(['cosine', 'pearson', 'spearman', 'kendall']):
    haxarr[i] = sns.heatmap(dist[m], cmap='RdBu_r', vmin=clim[0], vmax=clim[1],
                            ax=axarr[i],
                            xticklabels=ticklabels, yticklabels=ticklabels,
                            annot=True, annot_kws={'fontsize': 4},
                            square=True, linewidths=.5, cbar_kws={"shrink": .5})
    axarr[i].set_title(m)

fig1.suptitle(title_str)
plt.show()
# %% MAKE FIG 2 - CLUSTERMAP
sns.set_context('paper')
g = sns.clustermap(df_peakresp,
                   figsize=(8.5, 11), metric='cosine',
                   dendrogram_ratio=(.4, .2),
                   vmin=0.0, vmax=3.0,
                   cbar_pos=(.03, .2, .02, .2))

box = g.ax_heatmap.get_position()
cbox = g.ax_cbar.get_position()
g.ax_cbar.set_position([cbox.x0, box.y0, .02, .2])
plt.suptitle(f"date:{expt['date']}, fly={expt['fly']}, movie={expt['movie']}, conc={movie_conc}"
             f"\nscope={expt['microscope']}, reared=1-5ol+1-6ol",
             fontname='Myriad Pro')
fig2 = g.fig
plt.show()

# %% fig 3- ortho heatmap - for viewing cell clusters
importlib.reload(vis)
sns.set_context('paper')
#isort = np.argsort(embedding['embedding'].squeeze())
isort = proc['model'].isort
scopemap = df_frametimeinfo.scope_pulse_idx.to_numpy().reshape(1, -1)
olfmap = df_frametimeinfo.olf_pulse_idx.to_numpy().reshape(1, -1)
clustmap = isort[:, np.newaxis]

fig3, [ax, [ax_x_scope, ax_x_olf], ax_y] = vis.ortho_pcolormesh(dfof[isort, :], scope_pulse_idx=scopemap,
                                                                olf_pulse_idx=olfmap,
                                                                cell_clusters=clustmap)
fig3.suptitle(f"date={expt['date']}, fly={expt['fly']}, movie={expt['movie']}, conc={movie_conc}\n"
              f"scope={expt['microscope']}, reared=1-5ol+1-6ol, DATA=peakresp(filtdat['dfof'])",
              fontname='DejaVu Sans')
plt.show()
# %%
pdfname = PROC_DIR.joinpath('response_analysis(1).pdf')
with PdfPages(pdfname) as pdf:
    pdf.savefig(figure=fig0)
    pdf.savefig(figure=fig1)
    pdf.savefig(figure=fig2)
    pdf.savefig(figure=fig3)

    d = pdf.infodict()
    d['Title'] = 'Multipage PDF Example'
    d['Author'] = 'remy yang'
    d['Project'] = 'mb-odor-rearing'
    d['Data.date'] = expt['date']
    d['Data.fly'] = expt['fly']
    d['Data.movie'] = expt['movie']
    d['CreationDate'] = datetime.now()
    d['ModDate'] = datetime.today()
print('done saving to pdf.')

# %% plot rastermap w/ lines


# %%
fig, ax = vis.plot_rastermap(Fr, model.isort, title=f"Neural data embedding (dfof, z-scored)\n{expt['movie']}",
                             clim=(0, 1.0), box_aspect=1 / 3, stimframes=stim_frames, blockframes=block_frames,
                             xid=model.xid)
# %% plot cell responses, trial by trial (grouped by stimulus)
n_trials = ti['num_stim']

trialtrace_idx = df_frametimeinfo['trialtrace_idx'].to_numpy()
trial_traces = [None] * n_trials

for i in range(0, n_trials):
    dff = sp_sorted[:, trialtrace_idx == i + 1]
    trial_traces[i] = dff
# %%
sort_idx = df_stimulus['stimname'].argsort().tolist()
# sort_idx = list(range(0,n_trials))
stimname = df_stimulus['stimname'].tolist()

f, axarr = plt.subplots(1, n_trials, figsize=(11, 5), sharey=True)
for i in range(0, n_trials):
    im = axarr[i].imshow(trial_traces[sort_idx[i]], vmin=0.3, vmax=0.7, aspect=5, cmap='gray_r',
                         interpolation='none')
    axarr[i].set_title("\n".join(textwrap.wrap(stimname[sort_idx[i]], 8)), fontsize=7)

cbar = f.colorbar(axim, ax=axarr.ravel().tolist(), shrink=0.4)
f.suptitle('2021-04-23 fly 2 (1-5ol+1-6ol reared), movie_006 (-6)', fontsize=12)
plt.show()

# %%
f.savefig(PROC_DIR.joinpath('trials_Fc_rastersorted_odorgrouped.pdf'), bbox_inches="tight")

# %%
df_dat = pd.DataFrame(data=Fr)
df_stimcat = df_stimulus['stimname']
df_stimcat.index += 1
df_stimcat = pd.concat([pd.Series(['none']), df_stimcat])

mi = pd.MultiIndex.from_arrays([model.embedding.squeeze(), model.xid, np.arange(K) + 1],
                               names=('embedding', 'xid', 'ocid'))
df_dat = df_dat.set_index(mi, append=True)
df_dat.columns = pd.MultiIndex.from_frame(df_frametimeinfo.reset_index())
df_dat.head()
df_dat = df_dat.iloc[model.isort]
df_dat.head()

# %%
sns.set_style("darkgrid")
g = sns.JointGrid(height=4, ratio=20, space=0)
g.fig.set_figwidth(11)
g.fig.set_figheight(5)

# Plot by accessing the ax_joint, ax_marg_x, and ax_marg_y attributes

sns.heatmap(df_dat.values, ax=g.ax_joint, cmap='rocket', cbar=False, vmax=0.7, vmin=0.3,
            xticklabels=list(range(1, 1 + K)), yticklabels=False)

xid = df_dat.index.get_level_values('xid').to_numpy()  # .reshape(-1, 1)
lut = dict(zip(np.arange(8), sns.color_palette('Set2', 8)))
xid_colors = pd.Series(xid % 8).map(lut)

sns.heatmap(xid % 8, ax=g.ax_marg_y, vmin=0, vmax=50, cmap=sns.color_palette("Set2"), xticklabels=[None],
            yticklabels=[None],
            cbar=False)

df_fti = df_frametimeinfo.loc[:, ['olf_pulse_idx', 'scope_pulse_idx']].transpose()
sns.heatmap(df_fti, mask=df_fti.to_numpy() == 0,
            cmap=sns.hls_palette(16), vmin=1, vmax=16, ax=g.ax_marg_x, xticklabels=[None], cbar=False)
plt.show()

# %%
xid = df_dat.index.get_level_values('xid').to_numpy()  # .reshape(-1, 1)
lut = dict(zip(np.arange(8), sns.color_palette('Set2', 8)))
xid_colors = pd.Series(xid % 8).map(lut)

df_fti = df_frametimeinfo.loc[:, ['olf_pulse_idx', 'scope_pulse_idx']].transpose()
lut_time = dict(zip(np.arange(16), sns.hls_palette(16)))
time_colors = df_fti.applymap(lambda x: lut_time[x])
# %%

g = sns.clustermap(pd.DataFrame(df_dat.values), row_cluster=False, col_cluster=False,
                   row_colors=xid_colors.to_frame(), col_colors=time_colors,
                   linewidth=0, cmap='rocket', vmin=0.3, vmax=0.7, figsize=(11, 5),
                   colors_ratio=0.01, dendrogram_ratio=(0.1, 0.1))
plt.show()
# %%

# x = model.embedding[model.isort].flatten().tolist()
# clusters = fclusterdata(model.embedding[model.isort], t=20, criterion='maxclust')
clusters = model.xid[model.isort]

# network_pal = sns.cubehelix_palette(np.unique(clusters).size,
#                                     light=.9, dark=.1, reverse=True,
#                                     start=1, rot=-2)
# network_pal = sns.husl_palette(10)
# network_pal = sns.husl_palette(10, s=.45)
# sns.palplot(sns.hls_palette(10, h=.5))
network_pal = sns.hls_palette(100)

# n_unique = df_frametimeinfo['scope_pulse_idx'].unique().nonzero()[0].size
pal_list = ['Blues', 'seagreen']
name_list = ['scope_pulse_idx', 'olf_pulse_idx']
lut = {}
for nam, clr in zip(name_list, pal_list):
    uni = df_frametimeinfo[nam].unique().tolist()
    n_uni = len(uni)
    pal = sns.color_palette(clr, n_uni)
    lut = dict(zip(uni, pal))
sns.color_palette("Blues", np.unique(df_frametimeinfo['scope_pulse_idx'].to_numpy()).nonzero())
# temp = rng.shuffle(np.unique(clusters))
network_lut = dict(zip(arr, network_pal))
network_colors = pd.Series(clusters).map(network_lut)

g = sns.clustermap(dat, row_cluster=False, col_cluster=False, row_colors=network_colors, linewidth=0, xticklabels=False,
                   yticklabels=50, cmap='gray_r', vmin=0.3, vmax=0.7, figsize=[11, 5],
                   colors_ratio=0.01, dendrogram_ratio=(0.1, 0.1))
n_X = ops['n_X']
g.ax_heatmap.set_title(f'rastermap (n_X={n_X})')
g.cax.set_box_aspect(5)

linewidth = 1
linealpha = 0.75
if stim_frames is not None:
    [g.ax_heatmap.axvline(x=x, linewidth=linewidth, alpha=linealpha) for x in stim_frames]

if block_frames is not None:
    [g.ax_heatmap.axvline(x=x, color='k', linewidth=linewidth, alpha=linealpha) for x in block_frames]

plt.show()
# %%
grouped = dat.groupby(model.xid[model.isort], axis=0)
datmean = grouped.mean()
datmean.columns = pd.MultiIndex.from_frame(df_frametimeinfo)
datmean.head()
fig2 = plt.figure(figsize=(11, 5))
fig2.tight_layout()
ax = sns.heatmap(datmean, vmin=0.3, vmax=0.7, cmap='rocket', xticklabels=300, yticklabels=True,
                 cbar_kws={'orientation': 'vertical', 'shrink': 0.5})
ax.set_box_aspect(1 / 2)

linewidth = 1
linealpha = 0.5
if stim_frames is not None:
    [ax.axvline(x=x, color='w', linestyle=':', linewidth=linewidth, alpha=linealpha) for x in stim_frames]

if block_frames is not None:
    [ax.axvline(x=x, color='w', linewidth=linewidth, alpha=linealpha) for x in block_frames]

cbar = ax.collections[0].colorbar
plt.show()
# %%
trial = datmean.iloc[:, datmean.columns.get_level_values(6) == 1]

ax = sns.heatmap(trial, vmin=0.3, vmax=0.7, cmap='rocket', xticklabels=False, yticklabels=True,
                 cbar_kws={'orientation': 'vertical', 'shrink': 0.5})
ax.set_box_aspect(2)
plt.show()
# %%
fig2 = plt.figure()
df = pd.DataFrame.from_dict({'x': model.embedding.squeeze(),
                             'cluster': model.xid})
sns.rugplot(data=df, y="x", hue="cluster", palette='Set2')
plt.show()
# %%


plt.show()
# %%
widths = [20, 1]
heights = [1, 20]
gs_kw = dict(width_ratios=widths, height_ratios=heights)
fig6, f6_axes = plt.subplots(ncols=2, nrows=2,  # constrained_layout=True,
                             gridspec_kw=gs_kw)
f6_axes[1, 0].set_box_aspect(1 / 4)

ax = sns.heatmap(pd.DataFrame(np.arange(40)), ax=f6_axes[1, 1], vmin=0, vmax=50, cmap=network_pal, xticklabels=[None],
                 yticklabels=[None],
                 cbar=False, square=True)

sns.heatmap(df_fti, mask=df_fti.to_numpy() == 0,
            cmap=sns.hls_palette(16), vmin=1, vmax=16, ax=f6_axes[0, 0], xticklabels=[None], cbar=False)

sns.heatmap(datmean, ax=f6_axes[1, 0], vmin=0.3, vmax=0.7, cmap='rocket', xticklabels=[None], yticklabels=True,
            cbar=False)

plt.show()

# %%
X = pdist(model.embedding[model.isort])
Z = linkage(X, 'single')
dn = dendrogram(Z)

plt.show()

# %%
orn_responses = orns.orns()

f, ax = plt.subplots(figsize=(14, 10))
sns.heatmap(orn_responses, cmap='coolwarm', linewidths=0.05, ax=ax)
ax.tick_params(axis='y', labelsize=8)  # y-axis
ax.tick_params(axis='x', labelsize=8)  # y-axis
plt.show()

# %%
sns.set_style("whitegrid")
ax = sns.heatmap(orn_responses, vmin=0, vmax=288, linewidths=0.1, yticklabels=1)
plt.show()

plt.savefig('temp.eps')
