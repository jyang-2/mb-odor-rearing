from pathlib import Path
from collections import namedtuple
from dataclasses import dataclass, field
from typing import List
from pathlib import Path
import pandas as pd
import numpy as np
import suite2p
from suite2p.extraction import dcnv
from rearing import hallem, io
from rearing.ca import respvec, trace, vis
from scipy import stats, signal
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import svm
from sklearn import preprocessing
from sklearn.metrics.pairwise import cosine_similarity
import math
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from sklearn.svm import SVC
from scipy.signal import convolve
from dataload import SparsityAnalysis, Fly, Experiment, load_all_flydata, s2p_frametrim, adjust_for_bad_frames, stim2cat
import copy
from scipy.stats import zscore
from scipy.ndimage import filters
from scipy.cluster.hierarchy import dendrogram, linkage, fclusterdata
import importlib
import textwrap

# set plotting parameters
params = {'legend.fontsize': 8,
          'axes.labelsize': 8,
          'axes.titlesize': 10,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'font.family': ['sans-serif'],
          'figure.dpi': 200.0,
          'lines.linewidth': 1
          }
plt.rcParams.update(params)

Hallem = hallem.HallemSet()

PRJ_DIR = Path.cwd()
FIG_DIR = PRJ_DIR.joinpath('figures', 'fig_compare_response_strength')

db = pd.read_csv(Path.cwd().joinpath('src', 'db.csv'))
fly_list = load_all_flydata(db)


# %% select fly, experiment
fly = copy.deepcopy(fly_list[0])
mov = copy.deepcopy(fly.expt_list[0])

df_frametimeinfo0 = mov.df_frametimeinfo
df_stimulus0 = mov.df_stimulus

idx = (df_frametimeinfo0['isbadframe'] == 0).to_numpy()

df_frametimeinfo, df_stimulus, s2p_dat = adjust_for_bad_frames(df_frametimeinfo0, df_stimulus0, mov.s2p_dat)
df_stimulus = stim2cat(df_stimulus)

# %%
# ---------------------------------
# lowpass filter traces
# ---------------------------------
filt_ops = {'fs': s2p_dat['ops']['fs'],
            'cutoff': 0.5,
            'trans_width': 0.1,
            'numtaps': 160}

b = trace.get_lowpass_filter(filt_ops)

blocked_traces = trace.split_traces_to_blocks(s2p_dat['Fc'], df_frametimeinfo.scope_pulse_idx)
filt_blocked_traces = [signal.filtfilt(b, 1, y) for y in blocked_traces]
baseline = [trace.percentile_baseline(y, fs=s2p_dat['ops']['fs'], win=45, percentile=10) for y in filt_blocked_traces]
df = [y - b for y, b in zip(filt_blocked_traces, baseline)]
dfof0 = [y / b for y, b in zip(df, baseline)]
dfof = np.hstack(tuple(dfof0).astype('float32'))
#%%
# trace deconvolution

ops = mov.s2p_dat['ops']

# baseline operation, preprocess
Fc = dcnv.preprocess(
     F=s2p_dat['Fc'],
     baseline=ops['baseline'],
     win_baseline=ops['win_baseline'],
     sig_baseline=10,
     fs=ops['fs'],
     prctile_baseline=ops['prctile_baseline']
 )

# get spikes
spks = dcnv.oasis(F=Fc, batch_size=ops['batch_size'], tau=ops['tau'], fs=ops['fs'])

# get convolved calcium trace
efilt = np.exp(- np.linspace(0, 50, 200) / (ops['tau'] * ops['fs']))
efilt /= efilt.sum()
efilt = efilt[np.newaxis, :]
sout = convolve(spks, efilt, mode='same')

blocked_sout = trace.split_traces_to_blocks(sout, df_frametimeinfo.scope_pulse_idx)
#%%
win = int(3 * ops['fs'])
Fcs = filters.gaussian_filter1d(Fc, win)
#%%
cid=0
plt.plot(Fc[cid, :])
plt.plot(Fcs[cid, :], color='k')
plt.show()
#%%
plt.plot(sout)
plt.show()

plt.plot(Fc[0,:])
plt.show()
# %%
# compute responses, store in dict "responses"

df_stimulus_good = df_stimulus.loc[df_stimulus['isgoodtrial']]
stim_list = df_stimulus_good.stim_idx.to_numpy()

peakresp_nmean, baseline_avg, baseline_std = respvec.compute_response_vectors(dfof,
                                                                              stim_list=stim_list,
                                                                              baseline_idx=df_frametimeinfo.baseline_idx,
                                                                              peakresp_idx=df_frametimeinfo.peakresp_idx,
                                                                              maxk=3,
                                                                              stack_numpy=True)
peakresp_amp = peakresp_nmean - baseline_avg
responses = {'peakresp_mean': peakresp_nmean, 'baseline_avg': baseline_avg, 'baseline_std': baseline_std,
             'peakresp_amp': peakresp_amp}
#%%
df_peakresp = pd.DataFrame(peakresp_amp.transpose())
df_peakresp = df_peakresp.set_index(pd.MultiIndex.from_frame(df_stimulus_good.loc[:, ['stim_idx', 'stimname', 'conc']]))
df_peakresp.sort_index(level=['stimname', 'conc'], inplace=True)
df_peakresp = df_peakresp.transpose()



# %%

stimidx_list = df_stimulus_good['stim_idx'].to_numpy()
stimname_list = df_stimulus_good['stimname'].to_numpy()
sort_idx = df_stimulus_good['stimname'].sort_values().index.to_numpy()

df_peakresp = pd.DataFrame(peakresp_amp)
df_peakresp.columns = pd.MultiIndex.from_arrays([stimidx_list, stimname_list], names=('stim_idx', 'stimname'))
df_peakresp = df_peakresp.iloc[:, sort_idx]
#%%
dist = dict()
dist['pearson'] = df_peakresp.corr(method='pearson')
dist['spearman'] = df_peakresp.corr(method='spearman')
dist['kendall'] = df_peakresp.corr(method='kendall')
df_cosine_corr = pd.DataFrame(1 - squareform(pdist(df_peakresp.transpose(), metric='cosine')))
df_cosine_corr.index = dist['pearson'].index
df_cosine_corr.columns = dist['pearson'].columns
dist['cosine'] = df_cosine_corr
#%%        pdb.set_trace()

fig_1, axarr_1, haxarr_11 = vis.corrmats(dist, include_trialnum=False)
plt.show()
# %%
fig, axarr = plt.subplots(ncols=len(blocked_traces), nrows=1, constrained_layout=True, sharey=True, sharex=True,
                          figsize=(6, 6), dpi=200)
for y, ax in zip(blocked_dfof, axarr.flat):
    K, T = y.shape
    offset = 0
    for yy in y:
        ax.plot(yy + offset, linewidth=0.5)
        offset = offset - 1
    ax.set_ylim([-1 * K * 1, 1.5 * 1])
plt.show()
#%%
save_multipage_pdf(mov.PROC_DIR().joinpath('response_analysis.pdf'), [fig_1, fig])

#%%
cid = 5
blk = 1
fig, ax1 = plt.subplots()
ax1.plot(blocked_traces[blk][cid, :])
ax1.plot(filt_blocked_traces[blk][cid, :], color='k')
ax2 = ax1.twinx()
ax2.plot(blocked_dfof[blk][cid, :], alpha=0.7, color='r')
plt.show()
#%%
cid = 5
blk = 1
fig, ax1 = plt.subplots()
ax1.plot(s2p_dat['Fc'][cid, :], alpha=0.75)
ax1.plot(filt_traces[cid, :], color='k')
ax1.plot(np.concatenate(baseline, axis=1)[cid, :], color='k')

ax2 = ax1.twinx()
ax2.plot(dfof[cid, :], alpha=0.7, color='r')

plt.show()