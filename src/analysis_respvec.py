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
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore
from scipy.ndimage import filters
from scipy.cluster.hierarchy import dendrogram, linkage, fclusterdata
import importlib
import textwrap
from rastermap.mapping import Rastermap
from itertools import product
from rearing import hallem, io, meta
from rearing.ca import vis, trace, respvec
from sklearn import preprocessing
from datetime import datetime
import heapq
from functools import partialmethod
from collections import Counter


pp = pprint.PrettyPrinter(indent=2, compact=True)
Hallem = hallem.HallemSet()

params = {'legend.fontsize': 8,
          'axes.labelsize': 8,
          'axes.titlesize': 10,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'font.family': ['sans-serif'],
          'figure.dpi': 100.0
          }
plt.rcParams.update(params)

  



def load_s2p_files(fname_):
    vals = suite2p.gui.io.load_files(fname_)
    keys = ('stat', 'ops', 'Fcell', 'Fneu', 'Spks', 'iscell', 'probcell', 'redcell', 'probredcell', 'hasred')
    s2p_dat = dict(zip(keys, vals))
    return s2p_dat


def preprocess_s2p_dat(data_dict, neu_coeff=0.7):
    iscell = data_dict['iscell']
    F = data_dict['Fcell'][iscell, :]
    Fneu = data_dict['Fneu'][iscell, :]
    Fc = F - neu_coeff * Fneu
    dat = data_dict.copy()
    dat.update({'iscell': iscell, 'F': F, 'Fneu': Fneu, 'Fc': Fc})
    return dat


def grouped_mean_sparsity(sparsity_, index_):
    n_sparsity, n_trials = sparsity_.shape
    keys = np.unique(index_)
    d_ = dict()
    for k in keys:
        idx = np.argwhere(index_ == k)
        d_[k] = np.mean(sparsity_[:, idx], axis=1)
    return d_


def make_sparsity_long(sparsity_, n_std, stim_idx, stimname):
    df_sparsity = pd.DataFrame(sparsity_, columns=stim_idx)
    df_sparsity['n_std'] = n_std

    df_sparsity_long = df_sparsity.melt(id_vars='n_std', var_name='stim_idx', value_name='sparsity')

    d = dict(zip(stim_idx, stimname))
    df_sparsity_long['stimname'] = df_sparsity_long['stim_idx'].map(d)
    return df_sparsity_long


def make_fig0a(df, col_order=None):
    g = sns.FacetGrid(df, row="zeroed", col="stimname", hue='stimname', col_order=col_order,
                      aspect=0.8, height=4, sharex=True, sharey=False)
    g.map(sns.lineplot, "n_std", "sparsity")
    g.add_legend()
    g.map(plt.axhline, y=0.1, ls=":", c=".5")
    g.map(plt.axhline, y=0.05, ls=":", c=".5")

    for ax_ in g.axes[0]:
        ax_.set_ylim([-0.1, 1.0])
    for ax_ in g.axes[1]:
        ax_.set_ylim([-0.1, 0.4])
    g.fig.subplots_adjust(top=0.8)  # adjust the Figure in rp
    return g


def make_fig0b(df):
    g = sns.FacetGrid(df, col="zeroed", hue='stimname', aspect=0.8, height=4, sharey=False)
    g.map(sns.lineplot, "n_std", "sparsity")
    g.add_legend()
    g.map(plt.axhline, y=0.1, ls=":", c=".5")
    g.map(plt.axhline, y=0.05, ls=":", c=".5")
    g.axes[0, 0].set_ylim([-0.1, 1.0])
    g.axes[0, 1].set_ylim([-0.1, 0.4])
    g.fig.subplots_adjust(top=0.8)  # adjust the Figure in rp
    return g


def make_fig1(dist, ticklabels=None, include_trialnum=True):
    df_corr = dist['pearson']
    if ticklabels is None:
        if include_trialnum:
            ticklabels = [f"{item[1]} ({item[0]})" for item in df_corr.columns.to_list()]
        else:
            ticklabels = [item[1] for item in df_corr.columns.to_list()]

    clim = [-1.0, 1.0]
    fig1, axarr = plt.subplots(2, 2, figsize=(8, 8), constrained_layout=True)
    axarr = axarr.flatten()
    haxarr = [None] * 4

    for i, m in enumerate(['cosine', 'pearson', 'spearman', 'kendall']):
        haxarr[i] = sns.heatmap(dist[m].values, cmap='RdBu_r', vmin=clim[0], vmax=clim[1],
                                ax=axarr[i],
                                xticklabels=ticklabels, yticklabels=ticklabels,
                                annot=False, annot_kws={'fontsize': 4},
                                square=True, linewidths=.25, cbar_kws={"shrink": .3})
        axarr[i].set_title(m)
    return fig1, axarr, haxarr


def save_multipage_pdf(pdfname, figarr, expt_):
    with PdfPages(pdfname) as pdf:
        for fig in figarr:
            pdf.savefig(figure=fig)

        d = pdf.infodict()
        d['Title'] = 'Multipage PDF Example'
        d['Author'] = 'remy yang'
        d['Project'] = 'mb-odor-rearing'
        d['Data.date'] = expt_['date']
        d['Data.fly'] = expt_['fly']
        d['Data.movie'] = expt_['movie']
        d['CreationDate'] = datetime.now()
        d['ModDate'] = datetime.today()
    print('done saving to pdf.')


# %% Set dataset info and directory locations
def main():
    expt = dict()
    expt['date'] = '2021-06-26'
    expt['fly'] = 1
    expt['microscope'] = 'galvo-res'

    # make directories
    PRJ_DIR = Path.cwd()
    FLY_DIR = PRJ_DIR.joinpath('processed_data', expt['date'], str(expt['fly']))

    # load metadata yaml
    metadata = io.load_yaml(FLY_DIR.joinpath("metadata.yaml"))
    acquisitions = metadata['Acquisitions'][f"fly{expt['fly']}"]
    all_movies = [item['thorimage'] for item in acquisitions]

    #movie_list = ['movie_005', 'movie_006', 'movie_007', 'movie_008']
    movie_list = ['movie_substack_1_to_5']
    for movie_name in movie_list:
        # ---------------------------------------------------
        # set movie/acquisition number and make directories
        # ---------------------------------------------------
        # iterate through movie names in dataset
        expt['movie'] = movie_name
        iacq = [item['thorimage'] for item in acquisitions].index(expt['movie'])

        # make movie directories
        PROC_DIR = FLY_DIR.joinpath(expt['movie'])
        S2P_DIR = PROC_DIR.joinpath('suite2p', 'plane0')
        # --------------------------
        #  LOADING
        # --------------------------
        # load metadata pandas dfs
        df_sessioninfo = pd.read_excel(PROC_DIR.joinpath('metadata.xlsx'), sheet_name='sessioninfo')
        df_stimulus = pd.read_excel(PROC_DIR.joinpath('metadata.xlsx'), sheet_name='stimulus')
        df_stimulus['stimname'] = pd.Categorical(df_stimulus['stimname'], ordered=True,
                                                 categories=['empty', 'pfo', '1-5ol', '1-6ol', '1-5ol+1-6ol', 'EP'])


        # load ti
        ti = io.load_ti(PROC_DIR.joinpath('ti.mat'))

        # load frametimeinfo
        df_frametimeinfo = pd.read_csv(PROC_DIR.joinpath('frametimeinfo.csv'))
        fti = df_frametimeinfo.to_records()

        # load suite2p files
        s2p_dat0 = load_s2p_files(S2P_DIR.joinpath('stat.npy'))
        s2p_dat = preprocess_s2p_dat(s2p_dat0)

        # load rastermodel data
        # rastermodel = np.load(S2P_DIR.joinpath('rastermodel.npy'), allow_pickle=True).item()

        K, T = s2p_dat['Fc'].shape
        use_ti = False
        if use_ti:
            n_trials = ti['num_stim']
            n_blocks = ti['num_blocks']
        else:
            n_trials = np.unique(fti.olf_pulse_idx[fti.olf_pulse_idx>0]).size
            n_blocks = np.unique(fti.scope_pulse_idx[fti.scope_pulse_idx>0]).size

        movie_conc = df_stimulus.loc[:, ['conc_a', 'conc_b']].min().min()

        dataset = dict()
        dataset['metadata'] = metadata
        dataset['info'] = {'expt': expt, 'DIR': [FLY_DIR, PROC_DIR, S2P_DIR],
                           'df_stimulus': df_stimulus, 'df_frametimeinfo': df_frametimeinfo, 'ti': ti, 'fti': fti,
                           'K': K, 'T': T, 'n_blocks': n_blocks, 'movie_conc': movie_conc}
        # ------------------------------------
        # lowpass filter traces
        # ------------------------------------
        fs = s2p_dat['ops']['fs']  # Sample rate, Hz
        cutoff = 0.5  # Desired cutoff frequency, Hz
        trans_width = .1  # Width of transition from pass band to stop band, Hz
        numtaps = 203  # Size of the FIR filter.

        lowpass_args = {'fs': fs, 'cutoff': cutoff, 'trans_width': trans_width, 'numtaps': numtaps}

        # lowpass filter Fc, block by block
        b = signal.remez(numtaps, [0, cutoff, cutoff + trans_width, 0.5 * fs], [1, 0], Hz=fs)
        trial_traces = np.array_split(s2p_dat['Fc'], n_blocks, axis=1)

        filt_trial_traces = [signal.filtfilt(b, 1, y) for y in trial_traces]
        baseline = [trace.percentile_baseline(y, fs=6.1, win=45, percentile=20) for y in filt_trial_traces]
        df = [y - b for y, b in zip(filt_trial_traces, baseline)]
        dfof = [y / b for y, b in zip(df, baseline)]
        dfof = np.hstack(tuple(dfof)).astype('float32')

        prct_baseline_args = {'fs': 6.1, 'win': 45, 'percentile': 10}

        # # lowpass filter Fc, block by block
        # trial_traces = np.array_split(s2p_dat['Fc'], ti['num_blocks'], axis=1)
        # filt_trial_traces = [signal.filtfilt(b, 1, y) for y in trial_traces]

        # -----------------------------
        # calculate df/f
        # -----------------------------


        #dfof = [trace.extract_dfof(y, fs=fs, win=45, percentile=10) for y in filt_trial_traces]
        #dfof = np.concatenate(dfof_byblock, axis=1)
        dataset['traces'] = {'lowpass_args': lowpass_args, 'prct_baseline_args': prct_baseline_args, 'dfof': dfof}
        # ------------------------------------
        # response vectors
        # ------------------------------------
        stimidx_list = df_stimulus['stim_idx'].to_numpy()
        stimname_list = df_stimulus['stimname'].to_numpy()
        sort_idx = df_stimulus['stimname'].sort_values().index.to_numpy()

        stimconc_list = [item for item in [0, -6, -5, -4, -3] for i in range(3)]

        peakresp_nmean, baseline_avg, baseline_std = respvec.compute_response_vectors(dfof,
                                                                                      baseline_idx=fti.baseline_idx,
                                                                                      peakresp_idx=fti.peakresp_idx,
                                                                                      maxk=3,
                                                                                      stack_numpy=True)
        peakresp_amp = peakresp_nmean - baseline_avg

        n_std = np.arange(0, 20, 0.1)
        sparsity = respvec.compute_sparsity(n_std, peakresp_amp, baseline_std)

        keys = ['stimidx_list', 'stimname_list', 'sort_idx',
                'peakresp_nmean', 'baseline_avg', 'baseline_std', 'peakresp_amp',
                'n_std', 'sparsity']
        vals = [stimidx_list, stimname_list, sort_idx,
                peakresp_nmean, baseline_avg, baseline_std, peakresp_amp,
                n_std, sparsity]
        dataset['response'] = {k: v for k, v in zip(keys, vals)}
        # -----------------------------------
        # rastermap options
        # -----------------------------------
        run_rastermap = False
        if run_rastermap:
            ops = {'n_components': 1, 'n_X': 40, 'alpha': 1., 'K': 1.,
                   'nPC': 200, 'constraints': 2, 'annealing': True, 'init': 'pca',
                   'start_time': 0, 'end_time': -1}

            model = Rastermap(n_components=ops['n_components'], n_X=ops['n_X'], nPC=ops['nPC'],
                              init=ops['init'], alpha=ops['alpha'], K=ops['K'],
                              constraints=ops['constraints'], annealing=ops['annealing'])
            model.fit(dfof)
            proc = {'model': model, 'ops': ops, 'activity_mode': 'filtered dfof', 'sp': dfof}
            np.save(S2P_DIR.joinpath("rastermodel.npy"), proc)
            print(f"rastermodel.npy saved for movie={expt['movie']}")

        grpmean_sparsity = grouped_mean_sparsity(sparsity, stimname_list)

        if 'empty' in stimname_list:
            sparsity0 = sparsity - grpmean_sparsity['empty']
        elif 'pfo' in stimname_list:
            sparsity0 = sparsity - grpmean_sparsity['pfo']

        # sparsity dataframe, long form
        df_sparsity = make_sparsity_long(sparsity,
                                         n_std=n_std,
                                         stim_idx=stimidx_list,
                                         stimname=stimname_list)
        df_stimulus['stim_idx'].apply(lambda x: stimconc_list[x - 1])
        df_stimulus['stimconc'] = df_stimulus['stim_idx'].apply(lambda x: stimconc_list[x-1])

        df_sparsity['zeroed'] = False


        df_sparsity0 = make_sparsity_long(sparsity0,
                                          n_std=n_std,
                                          stim_idx=stimidx_list,
                                          stimname=stimname_list)

        df_sparsity0['zeroed'] = True

        title_str = f"\ndate={expt['date']}, fly={expt['fly']}, movie={expt['movie']}, conc={movie_conc}\n" \
                    f"scope={expt['microscope']}, reared=1-5ol+1-6ol, DATA=peakresp(filtdat['dfof'])"

        # ------------------------------------
        # FIG 0: sparsity vs n_std
        # ------------------------------------
        sns.set_style('darkgrid')

        run_sparsity_plots = True
        if run_sparsity_plots:
            col_order = ['pfo', '1-5ol', '1-6ol', '1-5ol+1-6ol']

            df_fig0 = pd.concat([df_sparsity, df_sparsity0], axis=0)
            g_0a = make_fig0a(df_fig0, col_order=col_order)
            g_0a.fig.suptitle(title_str)
            plt.show()

            # fig0b plotting
            g_0b = make_fig0b(df_fig0)
            g_0b.fig.suptitle(title_str)
            plt.show()

        # ------------------------------------
        # FIG 1: correlation matrices
        # ------------------------------------
        df_peakresp = pd.DataFrame(peakresp_amp)
        df_peakresp.columns = pd.MultiIndex.from_arrays([stimidx_list, stimname_list], names=('stim_idx', 'stimname'))
        df_peakresp = df_peakresp.iloc[:, sort_idx]

        dist = dict()
        dist['pearson'] = df_peakresp.corr(method='pearson')
        dist['spearman'] = df_peakresp.corr(method='spearman')
        dist['kendall'] = df_peakresp.corr(method='kendall')
        df_cosine_corr = pd.DataFrame(1 - squareform(pdist(df_peakresp.transpose(), metric='cosine')))
        df_cosine_corr.index = dist['pearson'].index
        df_cosine_corr.columns = dist['pearson'].columns
        dist['cosine'] = df_cosine_corr
        #        pdb.set_trace()

        fig_1, axarr_1, haxarr_11 = make_fig1(dist)
        plt.show()

        # keep in mind for pdist using scipy.spatial.distance
        # >>> import itertools
        # >>> list(combinations(range(3),2))
        # [(0, 1), (0, 2), (1, 2)]

        # ------------------------------------
        # FIG 2: rastermap clustered pcolormesh w/ marginal colormesh
        # ------------------------------------

        m = dict(zip(range(1, n_trials + 1), [stimname_ord.index(item) + 1 for item in stimname_list]))
        m[0] = 0
        # np.array(list(map(lambda x: m[x], fti.olf_pulse_idx))).reshape(1, -1),

        isort = rastermodel['model'].isort
        scopemap = df_frametimeinfo.scope_pulse_idx.to_numpy().reshape(1, -1)
        olfmap = df_frametimeinfo.olf_pulse_idx.to_numpy().reshape(1, -1)
        clustmap = rastermodel['model'].xid[isort]

        # axarr3 = [ax, [ax_x_scope, ax_x_olf], ax_y]
        fig3a, axarr3a = vis.ortho_pcolormesh(dfof[isort, :],
                                              scope_pulse_idx=scopemap,
                                              olf_pulse_idx=olfmap,
                                              cell_clusters=clustmap[:, np.newaxis])
        fig3a.suptitle(title_str)
        plt.show()

        df_dfof = pd.DataFrame(dfof[isort, :])
        df_dfof_meanclust = df_dfof.groupby(clustmap, axis=0, as_index=True).mean()
        fig3b, axarr3b = vis.ortho_pcolormesh(df_dfof_meanclust.values,
                                              scope_pulse_idx=scopemap,
                                              olf_pulse_idx=olfmap,
                                              cell_clusters=df_dfof_meanclust.index.to_numpy()[:, np.newaxis])
        fig3b.suptitle(title_str)
        plt.show()

        figlist = [g_0a.fig, g_0b.fig, fig_1, fig3a, fig3b]
        save_multipage_pdf(PROC_DIR.joinpath('response_analysis.pdf'), figlist, expt)
        print(f"\t{expt['movie']} plots saved.")

if __name__ == '__main__':
    main()


