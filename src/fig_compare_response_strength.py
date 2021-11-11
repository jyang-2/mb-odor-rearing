from typing import List
import itertools
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats, signal
from scipy.spatial.distance import pdist, squareform
from sklearn import preprocessing
from sklearn import svm

from dataload import SparsityAnalysis, load_all_flydata, make_fly_list
from rearing import hallem, io, corrmat
from rearing.ca import respvec, trace, vis

Hallem = hallem.HallemSet()

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams.update({'legend.fontsize': 8,
                     'axes.labelsize': 8,
                     'axes.titlesize': 10,
                     'xtick.labelsize': 8,
                     'ytick.labelsize': 8,
                     'font.family': ['sans-serif'],
                     'figure.dpi': 400.0
                     })




def legend_without_duplicate_labels(ax0):
    handles, labels = ax0.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax0.legend(*zip(*unique), loc='lower left', fontsize=8, frameon=True, ncol=3)


def save_multipage_pdf(pdfname, figarr):
    with PdfPages(pdfname) as pdf:
        for fig_ in figarr:
            pdf.savefig(figure=fig_)
    print('done saving to pdf.')


def calc_movie_conc(df_stimulus_):
    mconc = df_stimulus_['conc'].to_numpy()
    mconc = np.unique(mconc[mconc < 0])[0].tolist()
    return mconc


def postproc_df_stimulus(df_stimulus0, good_stim_idx=None, stimname_order=None):
    """
    Edits df_stimulus that is directly loaded from file
        - stimuli can be excluded (if they occur during a df_frametimeinfo bad frame)
        - 'stimname' can be converted from a string to an ordinal categorical var.

    :param df_stimulus0:    df_stimulus before processing, just loaded from file
    :type df_stimulus0:     pd.DataFrame
    :param good_stim_idx:   stimulus indices to include (stimuli not included will be dropped from DataFrame)
    :type good_stim_idx:    List, np.ndarray
    :param stimname_order:  ordering, use to convert df_stimulus['stimname'] into an ordinal categorical variable
    :type stimname_order:   List
    :return:    df_stimulus_new - selected stimuli only, 'stimname' categorical
    :rtype:     pd.DataFrame
    """

    df_stimulus_new = df_stimulus0.set_index('stim_idx')

    if good_stim_idx is not None:
        df_stimulus_new = df_stimulus_new.loc[good_stim_idx]

    if stimname_order is not None:
        df_stimulus_new['stimname'] = pd.Categorical(df_stimulus_new['stimname'], ordered=True,
                                                     categories=stimname_order)

    df_stimulus_new.reset_index(inplace=True)
    return df_stimulus_new





def triu_df_corrmat(df, value_name='corr'):
    keep = np.triu(np.ones(df.shape), k=1).astype('bool').reshape(df.size)
    temp = df.droplevel('stim_idx', axis=1).droplevel('stim_idx', axis=0)
    dfs = temp.stack()[keep].rename(value_name).to_frame()
    dfs.index.set_names(names=['stimname_a', 'stimname_b'], inplace=True)
    dfs.reset_index(inplace=True)
    dfs['stimname'] = list(zip(dfs['stimname_a'], dfs['stimname_b']))
    dfs_meansorted = dfs.loc[:, ['stimname', value_name]].groupby(['stimname']).mean().sort_values(by=value_name)
    return dfs, dfs_meansorted

    # if monotonic:
    #     ord = dfs_meansorted.index.to_list()
    # else:
    #     ord = dfs['stimname'].unique().tolist()


def process_experiment(mov, dfof_ops, all_stim_ord=['pfo', '1-5ol', '1-6ol', 'EP', '1-6ol+EP', 'benz'],
                       save_outputs=False):
    """
    parameter sets:
        - dfof_ops
        - filt_ops
        - stim_ord = ['pfo', '1-5ol', '1-6ol', 'EP', '1-6ol+EP']
    intermediate outputs to save:
        - baseline
        - df
        - dfof
        - spaan_std
        - spaan_thr
        - dist --> .npz (
        - df_peakresp
        - df_peakresp_avg


    :param mov:
    :type mov:
    :param dfof_ops:
    :type dfof_ops:
    :param all_stim_ord:
    :type all_stim_ord:
    :param save_outputs:
    :type save_outputs:
    :return:
    :rtype:
    """

    print(f"\n\nprocessing.........{mov}")

    ti = mov.ti
    fs = trace.compute_fs(mov.fti.frame_times)  # compute frame rate
    s2p_dat = mov.s2p_dat  # suite2p data (dict)
    df_frametimeinfo0 = mov.df_frametimeinfo  # info for each frame
    df_fti = df_frametimeinfo0.loc[df_frametimeinfo0['isbadframe'] == 0]

    # --------------------------------------------------------------
    # DF_STIMULUS - make 'stimname' a categorical variable
    # --------------------------------------------------------------
    stim_list = df_fti.trialtrace_idx[df_fti.trialtrace_idx > 0].unique()  # stimuli to include
    stim_ord = [item for item in all_stim_ord if item in mov.df_stimulus['stimname'].unique()]  # stimname ordering
    df_stimulus = postproc_df_stimulus(mov.df_stimulus.copy(), good_stim_idx=stim_list, stimname_order=stim_ord)

    # --------------------------------------------------------------
    # compute movie concentration
    # --------------------------------------------------------------
    movie_conc = calc_movie_conc(df_stimulus)

    # --------------------------------------------------------------
    # drop bad frames, if any
    # --------------------------------------------------------------
    Fc = s2p_dat['Fc'][:, df_frametimeinfo0.isbadframe == 0]
    blocked_traces = trace.split_traces_to_blocks(Fc, df_fti.scope_pulse_idx)

    # --------------------------------------------------------------
    # lowpass filter traces
    # --------------------------------------------------------------
    filt_ops = dict(fs=s2p_dat['ops']['fs'], cutoff=0.5, trans_width=0.05, numtaps=230)  # << PARAM DICT
    b = signal.remez(filt_ops['numtaps'], [0, filt_ops['cutoff'], filt_ops['cutoff'] + filt_ops['trans_width'],
                                           0.5 * filt_ops['fs']], [1, 0], Hz=filt_ops['fs'])
    filt_blocked_traces = [signal.filtfilt(b, 1, y) for y in blocked_traces]

    # --------------------------------------------------------------
    # compute dF/F
    # --------------------------------------------------------------
    baseline = [trace.percentile_baseline(y, fs=fs, win=dfof_ops['win'], percentile=dfof_ops['percentile']) for y in
                filt_blocked_traces]
    df = [y - b for y, b in zip(filt_blocked_traces, baseline)]
    dfof = [y / b for y, b in zip(df, baseline)]
    dfof_z = [stats.zscore(y, axis=1) for y in dfof]
    dfof = np.hstack(tuple(dfof)).astype('float32')
    K, T = dfof.shape

    # --------------------------------------------------------------
    # compute responses, store in dict "responses"
    # --------------------------------------------------------------
    peakresp_nmean, baseline_avg, baseline_std = respvec.compute_response_vectors(dfof,
                                                                                  stim_list=stim_list,
                                                                                  baseline_idx=df_fti.baseline_idx,
                                                                                  peakresp_idx=df_fti.peakresp_idx,
                                                                                  maxk=3,
                                                                                  stack_numpy=True)
    peakresp_amp = peakresp_nmean - baseline_avg
    # peakresp_amp[peakresp_amp < 0] = 0
    # ------------------------------------------
    # CORR MATS - calculate odor distances
    # ------------------------------------------
    dist, df_peakresp = compute_corrmats(peakresp_amp, df_stimulus=df_stimulus)
    df_peakresp_avg = df_peakresp.groupby(axis=1, level='stimname').mean()

    dist_triu = dict()
    dist_triu_asc = dict()
    for k in dist.keys():
        dfs, dfs_meansorted = triu_df_corrmat(dist[k], value_name=k)
        dist_triu[k] = dfs
        dist_triu_asc[k] = dfs_meansorted

    responses = dict(peakresp_mean=peakresp_nmean, baseline_avg=baseline_avg, baseline_std=baseline_std,
                     peakresp_amp=peakresp_amp, df_peakresp=df_peakresp, df_peakresp_avg=df_peakresp_avg)

    # compute sparsity info, add df_melted to list
    trial_dict = df_stimulus.loc[:, ['stimname', 'conc', 'block_idx']].to_dict(orient='list')

    # std dev thresholding - sparsity calculation
    stimname = df_stimulus['stimname'].to_numpy()
    n_std = np.arange(0, 20, 0.1)
    sparsity_std = respvec.compute_sparsity(n_std, peakresp_amp, baseline_std)
    sparsity_std_0 = sparsity_std - np.mean(sparsity_std[:, stimname == 'pfo'], axis=1)[:, np.newaxis]
    spaan = SparsityAnalysis(sparsity_std, n_std, stim_list.tolist(), trial_dict)
    spaan_zeroed = SparsityAnalysis(sparsity_std_0, n_std, stim_list.tolist(), trial_dict)

    # amplitude - sparsity calculation
    thr_list = np.arange(0, 5.0, 0.1)
    sparsity_thr = respvec.compute_sparsity_thr(thr_list, peakresp_amp, baseline_std, std_thr=3)
    sparsity_thr_0 = sparsity_thr - np.mean(sparsity_thr[:, stimname == 'pfo'], axis=1)[:, np.newaxis]
    spaan_thr = SparsityAnalysis(sparsity_thr, thr_list, stim_list.tolist(), trial_dict)
    spaan_thr_zeroed = SparsityAnalysis(sparsity_thr_0, thr_list, stim_list.tolist(), trial_dict)

    sparsity = dict(spaan_std=spaan, spaan_std_zeroed=spaan_zeroed,
                    spaan_thr=spaan_thr, spaan_thr_zeroed=spaan_thr_zeroed)

    # ------------------------------------------
    # SAVE INTERMEDIATE ANALYSIS RESULTS
    # ------------------------------------------
    subdirnam = f"__win{dfof_ops['win']}__prctile{dfof_ops['percentile']}"
    savesubdir = mov.PROC_DIR().joinpath(subdirnam)
    if not savesubdir.exists():
        savesubdir.mkdir()

    params = dict(filt_ops=filt_ops, dfof_ops=dfof_ops, stim_ord=stim_ord, movie_conc=movie_conc)
    with open(savesubdir.joinpath('params.yaml'), 'w') as outfile:
        yaml.dump(params, outfile, default_flow_style=False)

    np.save(savesubdir.joinpath('params.npy'), params)
    np.save(savesubdir.joinpath('filt_ops.npy'), filt_ops)
    np.save(savesubdir.joinpath('dfof_ops.npy'), dfof_ops)

    np.save(savesubdir.joinpath('dfof.npy'), dfof)
    np.save(savesubdir.joinpath('df.npy'), df)
    np.save(savesubdir.joinpath('baseline.npy'), baseline)

    np.save(savesubdir.joinpath('responses.npy'), responses)
    np.save(savesubdir.joinpath('sparsity.npy'), sparsity)

    np.save(savesubdir.joinpath('distances.npy'), dist)
    # np.save(savesubdir.joinpath('dist_triu.npy'), dist_triu)
    # np.save(savesubdir.joinpath('dist_triu_asc.npy'), dist_triu_asc)

    print(f"\t\tfiles saved.")





def make_respmat_tidy(df, name=None):
    df = df.stack(level=[0, 1])
    df = df.to_frame(name=name)
    return df


# def make_correlation_tidy():
#     dfof_ops = dict(win=30, percentile=20)
#
#     PRJ_DIR = Path.cwd()
#     #FIG_DIR = PRJ_DIR.joinpath('figures', 'fig_compare_response_strength')
#
#     # load fake 'database' file
#     db = pd.read_csv(Path.cwd().joinpath('src', 'db.csv'))
#     db = db.loc[2:]
#     fly_list = load_all_flydata(db)
#
#     dist_triu_lists = dict(cosine=[], pearson=[], spearman=[], kendall=[])
#
#     for fly in fly_list:
#         for mov in fly.expt_list:
#             subdirnam = f"__win{dfof_ops['win']}__prctile{dfof_ops['percentile']}"
#             loadsubdir = mov.PROC_DIR().joinpath(subdirnam)
#
#             params, dfof, sparsity, responses, dist, dist_triu, dist_triu_asc = load_processed_data(loadsubdir)
#
#             df_stimulus = postproc_df_stimulus(mov.df_stimulus.copy(), stimname_order=params['stim_ord'])
#
#             for k in dist_triu.keys():
#                 dist_triu[k]['conc'] = params['movie_conc']
#                 dist_triu[k]['fly_num'] = fly.fly_num
#                 dist_triu_lists[k].append(dist_triu[k])
#
#     df_triu_pearson = pd.concat(dist_triu_lists['pearson'], axis=0)
#
#
#     # Show the conditional means
#     sns.pointplot(x="value", y="measurement", hue="species",
#                   data=iris, dodge=.532, join=False, palette="dark",
#                   markers="d", scale=.75, ci=None)
#
#     # within fly, all concentrations
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
#     sns.pointplot(data=df_    df_binresp = respvec.get_binresponse(procdat['responses']['df_peakresp'],
#                                          procdat['responses']['baseline_std'],
#                                          n_std=3, thr=0.5)
#     df_tuning = df_binresp.groupby(by='stimname', axis=1).any()triu_pearson.loc[df_triu_pearson['fly_num'] == 1], x='stimname', y='pearson', hue='conc', ax=ax1,
#                   dodge=True, join=False)
#     # sns.stripplot(data=df_triu_pearson.loc[df_triu_pearson['fly_num']==1], x="stimname", y="pearson", hue='conc', ax=ax1,
#     #               jitter=False, alpha=0.75, dodge=True)
#     ax1.set_ylim((-0.2, 1.0))
#     ax1.set_xticklabels(ax.get_xticklabels(), rotation=90)
#     ax1.set_title(fly.metadata['Fly info'][1]['reared'])
#
#     # sns.stripplot(data=df_triu_pearson.loc[df_triu_pearson['fly_num'] == 2], x="stimname", y="pearson", hue='conc',
#     #               ax=ax2,
#     #               jitter=False, alpha=0.75, dodge=True)
#     sns.pointplot(data=df_triu_pearson.loc[df_triu_pearson['fly_num'] == 2], x='stimname', y='pearson', hue='conc', ax=ax2,
#                   dodge=True, join=False)
#     ax2.set_ylim((-0.2, 1.0))
#     ax2.set_xticklabels(ax.get_xticklabels(), rotation=90)
#     ax2.set_title(fly.metadata['Fly info'][2]['reared'])
#
#     sns.despine()
#     plt.show()
#     plt.savefig('all_conc_pearson_scatter_ptplt.pdf')
#
#     # within concentration, multiple flies
#     fig, ax1 = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
#     sns.stripplot(data=df_triu_pearson.loc[df_triu_pearson['conc']==-5], x="stimname", y="pearson", hue='fly_num', ax=ax1,
#                   jitter=False, alpha=0.75, dodge=True)
#     ax1.set_title('conc=-5 (pearson corr)')
#     ax1.set_ylim((-0.2, 1.0))
#     ax1.set_xticklabels(ax.get_xticklabels(), rotation=90)
#     plt.show()
#     fig.savefig('pearson_scatter_conc-5.pdf')
#
#     fig, ax1 = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
#     sns.stripplot(data=df_triu_pearson.loc[df_triu_pearson['conc'] == -4], x="stimname", y="pearson", hue='fly_num',
#                   ax=ax1,
#                   jitter=False, alpha=0.75, dodge=True)
#     ax1.set_title('conc=-4 (pearson corr)')
#     ax1.set_ylim((-0.2, 1.0))
#     ax1.set_xticklabels(ax.get_xticklabels(), rotation=90)
#     plt.show()
#     fig.savefig('pearson_scatter_conc-4.pdf')
#
#     fig, ax1 = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
#     sns.stripplot(data=df_triu_pearson.loc[df_triu_pearson['conc'] == -3], x="stimname", y="pearson", hue='fly_num',
#                   ax=ax1,
#                   jitter=False, alpha=0.75, dodge=True)
#     ax1.set_title('conc=-3 (pearson corr)')
#     ax1.set_ylim((-0.2, 1.0))
#     ax1.set_xticklabels(ax.get_xticklabels(), rotation=90)
#     plt.show()
#     fig.savefig('pearson_scatter_conc-3.pdf')
#
#     # within concentration, multiple flies =========================================
#     fig, ax1 = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
#     sns.pointplot(data=df_triu_pearson.loc[df_triu_pearson['conc'] == -5], x="stimname", y="pearson", hue='fly_num',
#                   ax=ax1, dodge=True, join=False)
#     ax1.set_title('conc=-5 (pearson corr)')
#     ax1.set_ylim((-0.2, 1.0))
#     ax1.set_xticklabels(ax.get_xticklabels(), rotation=90)
#     plt.show()
#     fig.savefig('pearson_scatter_ptplt_conc-5.pdf')
#     fig.savefig('pearson_scatter_ptplt_conc-5.png')
#
#     fig, ax1 = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
#     sns.pointplot(data=df_triu_pearson.loc[df_triu_pearson['conc'] == -4], x="stimname", y="pearson", hue='fly_num',
#                   ax=ax1, dodge=True, join=False)
#     ax1.set_title('conc=-4 (pearson corr)')
#     ax1.set_ylim((-0.2, 1.0))
#     ax1.set_xticklabels(ax.get_xticklabels(), rotation=90)
#     plt.show()
#     fig.savefig('pearson_scatter_ptplt_conc-4.pdf')
#     fig.savefig('pearson_scatter_ptplt_conc-4.png')
#
#     fig, ax1 = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
#     sns.pointplot(data=df_triu_pearson.loc[df_triu_pearson['conc'] == -3], x="stimname", y="pearson", hue='fly_num',
#                   ax=ax1, dodge=True, join=False)
#     ax1.set_title('conc=-3 (pearson corr)')
#     ax1.set_ylim((-0.2, 1.0))
#     ax1.set_xticklabels(ax.get_xticklabels(), rotation=90)
#     plt.show()
#     fig.savefig('pearson_scatter_ptplt_conc-3.pdf')
#     fig.savefig('pearson_scatter_ptplt_conc-3.png')
#
#         # for k in dist_triu_asc.keys():
#         #     dist_triu_asc[k]['conc'] = params['movie_conc']
#         #     dist_triu_asc[k]['fly_num'] = fly.fly_num
#         #
#
#         tstr = f"date: {fly.date}, fly: {fly.fly_num}, movies: {mov.movie} [{params['movie_conc']}] " \
#                f"\nrearing: {fly.metadata['Fly info'][fly.fly_num]['reared']}" \
#                f"\nwin={dfof_ops['win']}, prctile={dfof_ops['percentile']}"



def run_basic():
    db = io.load_db()
    fly_list = make_fly_list(db)

    fly = fly_list[2]
    mov = fly.expt_list[1]

    subdirlist = [item.name for item in mov.PROC_DIR().glob('__*')]
    loadsubdir = subdirlist[0]

    procdat = io.load_processed_data(mov.PROC_DIR().joinpath(loadsubdir))

    sparsity = procdat['sparsity']
    sparsity['spaan_std'].df_melted.to_csv(mov.PROC_DIR().joinpath(loadsubdir, 'df_spaan_std_tidy.csv'))


# %%

def main():
    PRJ_DIR = Path.cwd()
    FIG_DIR = PRJ_DIR.joinpath('figures', 'fig_compare_response_strength')

    # load fake 'database' file
    db = pd.read_csv(Path.cwd().joinpath('src', 'db.csv'))
    db = db.loc[2:]
    fly_list = load_all_flydata(db)

    # loop through flies, movies, parameters
    win_list = [45, 60, 30]
    prctile_list = [20, 45]
    # win_list = [30]
    # prctile_list = [20, 45]
    param_list = list(itertools.product(win_list, prctile_list))

    do_analysis = True
    do_datavis = False
    write_plots = False

    # dfof_ops = dict(win=30, percentile=20)
    for fly in fly_list:
        for mov in fly.expt_list:
            for item in param_list:
                dfof_ops = dict(win=item[0], percentile=item[1])

                if do_analysis:  # run analysis, save outputs
                    process_experiment(mov, dfof_ops)

                if do_datavis:
                    subdirnam = f"__win{dfof_ops['win']}__prctile{dfof_ops['percentile']}"
                    loadsubdir = mov.PROC_DIR().joinpath(subdirnam)
                    params, dfof, sparsity, responses, dist, dist_triu, dist_triu_asc = load_processed_data(loadsubdir)

                    df_stimulus = postproc_df_stimulus(mov.df_stimulus.copy(), stimname_order=params['stim_ord'])

                    tstr = f"date: {fly.date}, fly: {fly.fly_num}, movies: {mov.movie} [{params['movie_conc']}] " \
                           f"\nrearing: {fly.metadata['Fly info'][fly.fly_num]['reared']}" \
                           f"\nwin={dfof_ops['win']}, prctile={dfof_ops['percentile']}"

                    # dist_triu = dict()
                    # dist_triu_asc = dict()
                    # for k in dist.keys():
                    #     dfs, dfs_meansorted = triu_df_corrmat(dist[k], value_name=k)
                    #     dist_triu[k] = dfs
                    #     dist_triu_asc[k] = dfs_meansorted

                    fig_list = []

                    # CORR MATS
                    fig_corrmat, axarr_corrmat, haxarr_corrmat = vis.corrmats(dist)
                    fig_corrmat.suptitle(tstr)
                    fig_list.append(fig_corrmat)
                    plt.show()

                    # CORR SCATTER
                    # ----------------------------------
                    fig_corrscatter_pearson, ax = vis.corrscatter(dist['pearson'])
                    ax.set_title('pearson')
                    fig_corrscatter_pearson.suptitle(tstr)
                    fig_list.append(fig_corrscatter_pearson)
                    plt.show()

                    fig_corrscatter_cosine, ax = vis.corrscatter(dist['cosine'])
                    ax.set_title('cosine')
                    fig_corrscatter_cosine.suptitle(tstr)
                    fig_list.append(fig_corrscatter_cosine)
                    plt.show()

                    fig_corrscatter_spearman, ax = vis.corrscatter(dist['spearman'])
                    ax.set_title('spearman')
                    fig_corrscatter_spearman.suptitle(tstr)
                    fig_list.append(fig_corrscatter_spearman)
                    plt.show()

                    fig_corrscatter_kendall, ax = vis.corrscatter(dist['kendall'])
                    ax.set_title('kendall')
                    fig_corrscatter_kendall.suptitle(tstr)
                    fig_list.append(fig_corrscatter_kendall)
                    plt.show()

                    # PEAKRESP_AMP DISPLOT
                    # ----------------------------------
                    g = sns.displot(responses['df_peakresp_avg'].melt(value_name='peakresp_amp'), x="peakresp_amp",
                                    row='stimname', hue='stimname', kde=True, stat="density",
                                    height=2, aspect=3, facet_kws=dict(margin_titles=True))
                    fig_peakamp_dist = g.fig
                    fig_peakamp_dist.suptitle(tstr)
                    fig_list.append(fig_peakamp_dist)
                    plt.show()

                    # PAIRPLOT - pair plot of mean amplitudes
                    g = sns.pairplot(responses['df_peakresp_avg'], diag_kind="kde")
                    g.set(xlim=(-0.5, 5.5), ylim=(-0.5, 5.5))
                    g.fig.suptitle(tstr)
                    plt.show()
                    fig_peakamp_pairplot = g.fig
                    fig_list.append(fig_peakamp_pairplot)
                    plt.show()

                    # SPARSITY_STD
                    g = vis.make_sparsity_plot(sparsity['spaan_std'].df_melted, stim_order=params['stim_ord'],
                                               title=tstr)
                    fig_sparsity_std = g.fig
                    fig_list.append(fig_sparsity_std)

                    g = vis.make_sparsity_plot(sparsity['spaan_std_zeroed'].df_melted, stim_order=params['stim_ord'],
                                               title=tstr)
                    fig_sparsity_std0 = g.fig
                    fig_list.append(fig_sparsity_std0)
                    plt.show()

                    sns.set_style("ticks")
                    fig_sparsity_std_all, axarr = plt.subplots(1, 2, figsize=(8, 4))
                    sns.lineplot(data=sparsity['spaan_std'].df_melted, x='n_std', y='sparsity', hue='stimname',
                                 ax=axarr[0])
                    axarr[0].set_title('sparsity_std')
                    sns.lineplot(data=sparsity['spaan_std_zeroed'].df_melted, x='n_std', y='sparsity', hue='stimname',
                                 ax=axarr[1])
                    axarr[1].set_title('sparsity_std_zeroed')
                    fig_sparsity_std_all.suptitle(tstr)
                    sns.despine(offset=3)
                    fig_list.append(fig_sparsity_std_all)
                    plt.show()

                    # SPARSITY_THR
                    g = vis.make_sparsity_plot(sparsity['spaan_thr'].df_melted, stim_order=params['stim_ord'],
                                               title=tstr)
                    fig_sparsity_thr = g.fig
                    fig_list.append(fig_sparsity_thr)

                    g = vis.make_sparsity_plot(sparsity['spaan_thr_zeroed'].df_melted, stim_order=params['stim_ord'],
                                               title=tstr)
                    fig_sparsity_thr0 = g.fig
                    fig_list.append(fig_sparsity_thr0)
                    plt.show()

                    sns.set_style("ticks")
                    fig_sparsity_thr_all, axarr = plt.subplots(1, 2, figsize=(8, 4))
                    sns.lineplot(data=sparsity['spaan_thr'].df_melted, x='n_std', y='sparsity', hue='stimname',
                                 ax=axarr[0])
                    axarr[0].set_title('sparsity_thr')
                    axarr[0].set_xlabel('peakresp_amp')
                    sns.lineplot(data=sparsity['spaan_thr_zeroed'].df_melted, x='n_std', y='sparsity', hue='stimname',
                                 ax=axarr[1])
                    axarr[1].set_title('sparsity_thr_zeroed')
                    axarr[1].set_xlabel('peakresp_amp')
                    fig_sparsity_thr_all.suptitle(tstr)
                    sns.despine(offset=3)
                    fig_list.append(fig_sparsity_thr_all)
                    plt.show()

                    # write plots to multipage pdf file
                    if write_plots:
                        savenam = f"sparsity_analysis_win{dfof_ops['win']}_prctile{dfof_ops['percentile']}.pdf"
                        pdfname = loadsubdir.joinpath(savenam)
                        with PdfPages(pdfname) as pdf:
                            for f in fig_list:
                                pdf.savefig(figure=f)
                        print('done saving to pdf.')

                    plt.close('all')


# %%
# def make_plots():


# for mov in [fly.expt_list[0]]:
#     tstr = f"date: {fly.date}, fly: {fly.fly_num}, movies: {mov.movie} [{movie_conc}] " \
#            f"\nrearing: {fly.metadata['Fly info'][fly.fly_num]['reared']}" \
#            f"\nwin={dfof_ops['win']}, prctile={dfof_ops['percentile']}"

# fig_corrmat, axarr_corrmat, haxarr_corrmat = vis.corrmats(dist)
# fig_corrmat.suptitle(tstr)
# plt.show()

# # SAVE INTERMEDIATE ANALYSIS RESULTS
# subdirnam = f"win{dfof_ops['win']}__prctile{dfof_ops['percentile']}"
#
# savesubdir = mov.PROC_DIR().joinpath(subdirnam)
# if not savesubdir.exists():py
#     savesubdir.mkdir()
#
# save_intermediate_outputs = True
# if save_intermediate_outputs:
#     params = dict(filt_ops=filt_ops, dfof_ops=dfof_ops, stim_ord=stim_ord, movie_conc=movie_conc)
#     np.save(savesubdir.joinpath('params'), params)
#     np.save(savesubdir.joinpath('filt_ops'), filt_ops)
#     np.save(savesubdir.joinpath('dfof_ops.npy'), dfof_ops)
#     np.save(savesubdir.joinpath('dfof.npy'), dfof)
#     np.save(savesubdir.joinpath('responses.npy'), responses)
#     np.save(savesubdir.joinpath('sparsity.npy'), sparsity)

#
# # SPARSITY PLOTS
# # ----------------------------
# # sparsity plot a
# # ----------------------------
# # df_sparsity_all = pd.concat(df_list, axis=0)
# # stim_order = [item for item in ord if item in df_sparsity['stimname'].unique()]
# # stim_ord = ['pfo', '1-5ol', '1-6ol', 'EP', '1-6ol+EP']
# # g = vis.make_sparsity_plot(df_sparsity, stim_order=stim_ord, title=tstr)
# df_sparsity = spaan.df_melted
# g = vis.make_sparsity_plot(df_sparsity, stim_order=stim_ord, title=tstr)
# fig_sparsity_ai = g.fig
# plt.show()
#
# g = vis.make_sparsity_plot(spaan_zeroed.df_melted, stim_order=stim_ord, title=tstr)
# fig_sparsity_aii = g.fig
# plt.show()
#
# # ----------------------------
# # sparsity plot b
# # ----------------------------
# sns.set_style("ticks")
# fig_sparsity_b, ax = plt.subplots(figsize=(4, 4))
# sns.lineplot(data=df_sparsity, x='n_std', y='sparsity', hue='stimname', ax=ax)
# ax.set_title(tstr)
# sns.despine(offset=3)
# plt.show()
#
# # --------------------------------------
# # sparsity df/f amplitude threshold plot
# # --------------------------------------
# sns.set_style("ticks")
# fig_sparsity_ci, ax = plt.subplots(figsize=(4, 4))
# sns.lineplot(data=spaan_thr.df_melted, x='n_std', y='sparsity', hue='stimname', ax=ax)
# ax.set_title(tstr)
# ax.set_xlabel('peakamp threshold')
# ax.set_ylabel("response sparsity (fraction cells responding)")
# sns.despine(offset=3)
# plt.show()
#
# fig_sparsity_cii, ax = plt.subplots(figsize=(4, 4))
# sns.lineplot(data=spaan_thr_zeroed.df_melted, x='n_std', y='sparsity', hue='stimname', ax=ax)
# ax.set_title(tstr)
# ax.set_xlabel('peakamp threshold')
# ax.set_ylabel("response sparsity (fraction cells responding)")
# sns.despine(offset=3)
# plt.show()
#
# # pair plot of mean amplitudes
# g = sns.pairplot(df_peakresp_avg, diag_kind="kde")
# g.set(xlim=(-0.5, 5.5), ylim=(-0.5, 5.5))
# g.fig.suptitle(tstr)
# plt.show()
# fig_peakamp_pairplot = g.fig
#
# savenam = f"sparsity_analysis_win{dfof_ops['win']}_prctile{dfof_ops['percentile']}.pdf"
# savedir = mov.PROC_DIR().joinpath('traces')
# if not savedir.exists():
#     savedir.mkdir()
#     np.save(savedir.joinpath('filt_ops'), filt_ops)
#     np.save(savedir.joinpath('dfof.npy'), dfof)
#     np.save(savedir.joinpath('responses.npy'), responses)
#     np.save(savedir.joinpath('sparsity.npy'), sparsity)
#
# pdfname = savedir.joinpath(savenam)
# with PdfPages(pdfname) as pdf:
#     pdf.savefig(figure=fig_corrmat)
#     pdf.savefig(figure=fig_sparsity_ai)
#     pdf.savefig(figure=fig_sparsity_aii)
#     pdf.savefig(figure=fig_sparsity_b)
#     pdf.savefig(figure=fig_sparsity_ci)
#     pdf.savefig(figure=fig_sparsity_cii)
#     pdf.savefig(figure=fig_peakamp_pairplot)
# print('done saving to pdf.')
# %% response rate plot
respOpts = dict(thr=1.0, std_thr=3)
responders = (peakresp_amp > (respOpts['std_thr'] * baseline_std)) & (peakresp_amp > respOpts['thr'])
response_rate = np.sum(responders, axis=0) / responders.shape[0]
df_response_rate = df_stimulus.copy(deep=True)
df_response_rate['response_rate'] = response_rate.transpose()
df_response_rate.sort_values(by='stimname', inplace=True)

# Draw a nested barplot by species and sex
g = sns.catplot(
    data=df_response_rate, kind="bar",
    x="stimname", y="response_rate", hue="stimname",
    ci="sd")
g.despine()
g.set_axis_labels("odors", "response rate (sparsity)")

# %% SVM ANALYSIS, PLOTTING

# stim_frames: stim_on frame #

if fly.date == '2021-03-23' and fly.fly_num == 1:
    stim_fct = ti['stim_fct'][:15]
    stim_fct = stim_fct[1:]
else:
    stim_fct = ti['stim_fct']
stim_frames = np.interp(stim_fct, df_fti.frame_times, df_fti.index)
# stim_frames = stim_frames - stim_frames[0]
# %%
# df: dataframe w/ column = cell, row = stimulus, index = (odors, concentration)
# PREPROCESSING

df_stimulus = mov.df_stimulus.copy()
df_stimulus = df_stimulus.set_index('stim_idx').loc[stim_list]

training_type = 'timeseries'
preproc_type = 'none'

scaler = preprocessing.StandardScaler()
if training_type == 'respvec':
    if preproc_type == 'StandardScaler':
        data = scaler.fit_transform(peakresp_nmean.transpose())
    elif preproc_type == 'none':
        data = peakresp_nmean.transpose()
    df = pd.DataFrame(data,
                      index=list(zip(df_stimulus.stimname, df_stimulus.conc)))
elif training_type == 'timeseries':
    if preproc_type == 'StandardScaler':
        scaled_dfof = scaler.fit_transform(dfof.transpose())
        scope_pulse_idx = df_fti['train_pulse_idx'].to_numpy()
        data = scaled_dfof[scope_pulse_idx > 0, :]
    elif preproc_type == 'none':
        scaled_dfof = dfof.transpose()
        scope_pulse_idx = df_fti['train_pulse_idx'].to_numpy()
        data = scaled_dfof[scope_pulse_idx > 0, :]
    df_temp = df_stimulus.loc[:, ['stimname', 'conc']].loc[scope_pulse_idx[scope_pulse_idx > 0]]
    df = pd.DataFrame(data, index=list(zip(df_temp.stimname, df_temp.conc)))

# conc_list: list of concentrations to use for training
conc_list = [-6, -5, -4, -3]

# order stimnames
df_stimulus = mov.df_stimulus.copy()
df_stimulus = df_stimulus.set_index('stim_idx').loc[stim_list]

ord_all = ['empty', 'pfo', '1-5ol', '1-6ol', '1-5ol+1-6ol', 'EP']
ord_mov = list(filter(lambda x: x in df_stimulus.stimname.unique(), ord_all))
cat_type = pd.api.types.CategoricalDtype(categories=ord_mov, ordered=True)
df_stimulus['stimname'] = df_stimulus['stimname'].astype(cat_type)

stimnames = df_stimulus['stimname'].to_list()
concs = df_stimulus['conc'].to_list()

# %% PLOTTING - iterate through concentrations
fig_list = [None] * len(conc_list)
tstr = f"date: {fly.date}, fly: {fly.fly_num}, movies: {mov.movie}"

for ifig, svc_conc in enumerate(conc_list):
    # format data for sklearn.svm
    if (('1-5ol', svc_conc) in df.index.to_list()) and (('1-5ol', svc_conc) in df.index.to_list()):
        label_dict = {('1-5ol', svc_conc): 0, ('1-6ol', svc_conc): 1}
        X = df.loc[[('1-5ol', svc_conc), ('1-6ol', svc_conc)], :].to_numpy()  # [# stimuli] x [# cells]
        Y = df.loc[[('1-5ol', svc_conc), ('1-6ol', svc_conc)], :].index.map(label_dict).to_numpy()

        # -------------------------------------
        # fit svm
        # -------------------------------------
        svc = svm.SVC(kernel='linear')
        svc.fit(X, Y)

        if training_type == 'respvec':
            if preproc_type == 'none':
                dfof_prj = svc.decision_function(dfof.transpose())
            elif preproc_type == 'StandardScaler':
                dfof_prj = svc.decision_function(scaler.transform(dfof.transpose()))
            resp_prj = svc.decision_function(data)
        elif training_type == 'timeseries':
            dfof_prj = svc.decision_function(scaled_dfof)
            resp_prj = svc.decision_function(df.groupby(by=df_temp.index).mean().to_numpy())

        # -------------------------------------
        # SVM PROJECTION FIGURE
        # -------------------------------------
        fig_list[ifig], ax = plt.subplots(1, 2, figsize=(8, 4), sharey='all')
        # -------------------------------------
        # subplot 1, plot of dfof_prj w/ stim_on times indicated
        # -------------------------------------
        ax[0].plot(df_fti.index.to_numpy(), dfof_prj, label='decision mag.')

        cat = df_stimulus['stimname'].astype('category').cat.codes.to_list()

        # make colormap, 1 colormap per concentration
        cmap = sns.color_palette()
        clr = [cmap[item] for item in cat]

        # plot lines indicating start of stimulus, & stimname
        for x, cc, c, s in zip(stim_frames, concs, clr, stimnames):
            if svc_conc == cc:
                ax[0].axvline(x=x, c=c, zorder=1, alpha=1, label=s)
            else:
                ax[0].axvline(x=x, c=c, zorder=1, alpha=0.3, label=s)
            ax[0].text(x=x, y=0.5, s=f"{s} [{cc}]", rotation=90,
                       horizontalalignment='left', verticalalignment='bottom',
                       rotation_mode='anchor', fontsize=6, color='k')

        legend_without_duplicate_labels(ax[0])
        ax[0].set_xlabel('frame')
        ax[0].set_ylabel('svm decision')
        ax[0].set_title(f'trained on 1-5ol vs. 1-6ol at conc={svc_conc}')
        ax[0].set_ylim([-2, 2])

        # plot svm projection (distance to decision boundary)
        df_mag = pd.DataFrame(resp_prj, columns={'svm decision'})
        df_mag['stimname'] = stimnames
        df_mag['conc'] = concs
        g = sns.barplot(data=df_mag, x='stimname', y='svm decision', hue='conc', ax=ax[1],
                        order=ord,
                        palette=sns.color_palette())
        g.set_xticklabels(g.get_xticklabels(), rotation=45)
        movie_conc = df_stimulus.conc.loc[df_stimulus.conc < 0].unique()
        ttstr = f"SVM magnitude projection (movie={mov.movie}, conc = {movie_conc})\n" \
                + tstr + "\n" + f"train. type={training_type}, preproc. type={preproc_type}"
        fig_list[ifig].subplots_adjust(top=0.8, hspace=0.3, bottom=0.2)
        fig_list[ifig].suptitle(ttstr, fontsize=10)
        plt.show()

# %%


# %%
for ifig, svc_conc in enumerate(conc_list):
    # format data for sklearn.svm
    if (('1-5ol', svc_conc) in df.index.to_list()) and (('1-5ol', svc_conc) in df.index.to_list()):
        label_dict = {('1-5ol', svc_conc): 0, ('1-6ol', svc_conc): 1}
        # X1 =

        X = df.loc[[('1-5ol', svc_conc), ('1-6ol', svc_conc)], :].to_numpy()  # [# stimuli] x [# cells]
        Y = df.loc[[('1-5ol', svc_conc), ('1-6ol', svc_conc)], :].index.map(label_dict).to_numpy()

# %%
for ifig, svc_conc in enumerate(conc_list):
    # format data for sklearn.svm
    if (('1-5ol', svc_conc) in df.index.to_list()) and (('1-5ol', svc_conc) in df.index.to_list()):
        label_dict = {('1-5ol', svc_conc): 0, ('1-6ol', svc_conc): 1}
        X = df.loc[[('1-5ol', svc_conc), ('1-6ol', svc_conc)], :].to_numpy()  # [# stimuli] x [# cells]
        Y = df.loc[[('1-5ol', svc_conc), ('1-6ol', svc_conc)], :].index.map(label_dict).to_numpy()

        # -------------------------------------
        # fit svm
        # -------------------------------------
        svc = svm.SVC(kernel='linear')
        svc.fit(X, Y)

        if training_type == 'respvec':
            if preproc_type == 'none':
                dfof_prj = svc.decision_function(dfof.transpose())
            elif preproc_type == 'StandardScaler':
                dfof_prj = svc.decision_function(scaler.transform(dfof.transpose()))
            resp_prj = svc.decision_function(data)
        elif training_type == 'timeseries':
            dfof_prj = svc.decision_function(scaled_dfof)
            resp_prj = svc.decision_function(df.groupby(by=df_temp.index).mean().to_numpy())

        # -------------------------------------
        # SVM PROJECTION FIGURE
        # -------------------------------------
        fig_list[ifig], ax = plt.subplots(1, 2, figsize=(8, 4), sharey='all')
        # -------------------------------------
        # subplot 1, plot of dfof_prj w/ stim_on times indicated
        # -------------------------------------
        ax[0].plot(df_fti.index.to_numpy(), dfof_prj, label='decision mag.')

        cat = df_stimulus['stimname'].astype('category').cat.codes.to_list()

        # make colormap, 1 colormap per concentration
        cmap = sns.color_palette()
        clr = [cmap[item] for item in cat]

        # plot lines indicating start of stimulus, & stimname
        for x, cc, c, s in zip(stim_frames, concs, clr, stimnames):
            if svc_conc == cc:
                ax[0].axvline(x=x, c=c, zorder=1, alpha=1, label=s)
            else:
                ax[0].axvline(x=x, c=c, zorder=1, alpha=0.3, label=s)
            ax[0].text(x=x, y=0.5, s=f"{s} [{cc}]", rotation=90,
                       horizontalalignment='left', verticalalignment='bottom',
                       rotation_mode='anchor', fontsize=6, color='k')

        legend_without_duplicate_labels(ax[0])
        ax[0].set_xlabel('frame')
        ax[0].set_ylabel('svm decision')
        ax[0].set_title(f'trained on 1-5ol vs. 1-6ol at conc={svc_conc}')
        ax[0].set_ylim([-2, 2])

        # plot svm projection (distance to decision boundary)
        df_mag = pd.DataFrame(resp_prj, columns={'svm decision'})
        df_mag['stimname'] = stimnames
        df_mag['conc'] = concs
        g = sns.barplot(data=df_mag, x='stimname', y='svm decision', hue='conc', ax=ax[1],
                        order=ord,
                        palette=sns.color_palette())
        g.set_xticklabels(g.get_xticklabels(), rotation=45)
        movie_conc = df_stimulus.conc.loc[df_stimulus.conc < 0].unique()
        ttstr = f"SVM magnitude projection (movie={mov.movie}, conc = {movie_conc})\n" \
                + tstr + "\n" + f"train. type={training_type}, preproc. type={preproc_type}"
        fig_list[ifig].subplots_adjust(top=0.8, hspace=0.3, bottom=0.2)
        fig_list[ifig].suptitle(ttstr, fontsize=10)
        plt.show()

# %%
save_multipage_pdf(fly.FLY_DIR.joinpath('svm_projection_AvsB_StandardScaler_dfofTrained.pdf'), fig_list)

# %% sparsity plot

ord = ['pfo', '1-5ol', '1-6ol', '1-5ol+1-6ol']

# save to pdf
g.fig.savefig(fly.FLY_DIR.joinpath(f"sparsity__{fly.date}__{fly.fly_num}.pdf"))
#
