#  STEP 1. Setup

from pathlib import PureWindowsPath, PurePath, Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import copy
import json
from collections import namedtuple
from functools import partial

from scipy import signal
from scipy.spatial.distance import squareform, pdist
from scipy.stats import zscore
from scipy.ndimage import filters, gaussian_filter1d
from scipy.io import savemat

from ryim import dirutils, analysis

from rearing import tidy
# from rearing.ca import trace, vis, represp

import pyarrow.feather as feather

import plotly.express as px
import plotly.io as pio

pio.renderers.default = "browser"
import pprint as pp
from step00_process_suite2p import postprocess_metadata, load_metadata, get_title_str

# %%

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
#     _session_info = {key: value for key, value in _nt_meta_data.session_info.to_dict().items() if key.isidentifier()}
#
#     return MetaData(_df_stimulus,
#                     _df_frametimeinfo,
#                     _nt_meta_data.timestamps,
#                     _session_info,
#                     _nt_meta_data.xml_meta,
#                     _nt_meta_data.sync_meta,
#                     _nt_meta_data.meta_yaml
#                     )

def lowpass_filter(_cells_x_time, Fs, **kwargs):
    """ Lowpass filter fluorescence. Default cutoff freq. is 0.5 Hz.
    :param _cells_x_time:
    :param float fs: frame rate
    :param cutoff: default = 0.05
    :param trans_width: default = 0.05
    :param numtaps: default = 230

    :return:
        np.array : lowpass filtered fluorescence traces, cells_x_time
    """
    filt_ops = dict(cutoff=0.5, trans_width=0.05, numtaps=230)
    filt_ops.update(kwargs)

    b = signal.remez(filt_ops['numtaps'],
                     [0, filt_ops['cutoff'],
                      filt_ops['cutoff'] + filt_ops['trans_width'],
                      0.5 * Fs], [1, 0], Hz=Fs)
    return signal.filtfilt(b, 1, _cells_x_time)


def calculate_noise_levels(neurons_x_time, frame_rate):
    """"
    dNMF function, copied off github

    Computes the noise levels for each neuron of the input matrix 'dF_traces'.
    The noise level is computed as the median absolute dF/F difference
    between two subsequent time points. This is a outlier-robust measurement
    that converges to the simple standard deviation of the dF/F trace for
    uncorrelated and outlier-free dF/F traces.
    Afterwards, the value is divided by the square root of the frame rate
    in order to make it comparable across recordings with different frame rates.
    input: dF_traces (matrix with nb_neurons x time_points)
    output: vector of noise levels for all neurons
    """
    dF_traces = neurons_x_time
    noise_levels = np.nanmedian(np.abs(np.diff(dF_traces, axis=-1)), axis=-1) / np.sqrt(frame_rate)
    return noise_levels * 100  # scale noise levels to percent


def compute_corrmats(peakresp_amp_, df_stimulus=None):
    df_peakresp = pd.DataFrame(peakresp_amp_)

    if df_stimulus is not None:
        sort_idx = df_stimulus['stimname'].sort_values().index.to_list()
        df_peakresp.columns = pd.MultiIndex.from_frame(df_stimulus.loc[:, ['stim_idx', 'stimname']])
        df_peakresp = df_peakresp.iloc[:, sort_idx]

    dist = dict()
    dist['pearson'] = df_peakresp.corr(method='pearson')
    dist['spearman'] = df_peakresp.corr(method='spearman')
    dist['kendall'] = df_peakresp.corr(method='kendall')

    df_cosine_corr = pd.DataFrame(1 - squareform(pdist(df_peakresp.transpose(), metric='cosine')))
    df_cosine_corr.index = dist['pearson'].index
    df_cosine_corr.columns = dist['pearson'].columns
    dist['cosine'] = df_cosine_corr

    for k, df in dist.items():
        df.columns.names = ['stim_idx_col', 'stimname_col']
        df.index.names = ['stim_idx_row', 'stimname_row']

    return dist, df_peakresp


def mean_corrmat(df_corrmat):
    """
    Calculates grouped means for odor pairs in corrmat.

    Parameters
    ----------
    df_corrmat - correlation matrix in "dist" dictionary returned by respvec.compute_corrmats(...)

    Returns
    -------
        df_mean_corrman - pd.DataFrame, with index = stimname_row, and columns = 'stimname_col'

    """
    # wide to long data format
    df = tidy.tidy_df_corrmat(df_corrmat)

    # set diagonals to np.nan to ignore odor identity distances (correlation between an odor trial and itself)
    df.loc[df['stim_idx_row'] == df['stim_idx_col']] = np.nan

    # calculate grouped means of odor pairs, and pivot to wide format,
    # where df_mean_corrmat.values[i, j] = correlation(trial i, trial j)
    df_mean_corrmat = df.loc[:, ['stimname_row', 'stimname_col', 'value']] \
        .groupby(by=['stimname_row', 'stimname_col']) \
        .mean() \
        .pivot_table(index='stimname_row', columns='stimname_col', values='value')

    return df_mean_corrmat


def mean_corrmats(dist):
    """
    Compute odor pair correlation means, for ever corrmat in dist

    Parameters
    ----------
    dist : dictionary w/ different corrmats returned by respvec.compute_corrmats(...)

    Returns
    -------
        mean_dist : dict w/ keys {'cosine', 'pearson', 'kendall', 'spearman'}

    """
    mean_dist = dict()
    for k, v in dist.items():
        mean_dist[k] = mean_corrmat(v)
    return mean_dist


def reflect_tidy_corr(tidy_df):
    """
    Make tidy odor-odor pair correlation (long format pd.DataFrame) for facet violin plots

    Parameters
    ----------
    tidy_df

    Returns
    -------

    """
    tidy_copy = copy.deepcopy(tidy_df)
    tidy_copy['stimname_row'] = tidy_df['stimname_col']
    tidy_copy['stimname_col'] = tidy_df['stimname_row']
    tidy_copy['stim_idx_row'] = tidy_df['stim_idx_col']
    tidy_copy['stim_idx_col'] = tidy_df['stim_idx_row']
    tidy_copy = tidy_copy.loc[tidy_copy['stimname_row'] != tidy_copy['stimname_col'], :]

    return tidy_df.append(tidy_copy)


def faceted_corrscatter(f, metric):
    """

    Parameters
    ----------
    f : file in the same folder as
                - tidy_corr_all.feather and
                - tidy_corr_resp.feather

    metric : which correlation/pdist metric to use {'cosine', 'pearson', 'kendall', 'tau'}

    Returns
    -------
        g : sns.catplot, or a facetgrid of odor-odor correlation distances (trial by trial)
            to save the figure, use `g.fig.savefig('corrscatter.png')`

    """
    tidy_corr_all = feather.read_feather(f.with_name("tidy_corr_all.feather"))
    tidy_corr_resp = feather.read_feather(f.with_name("tidy_corr_resp.feather"))

    tidy_violin_all = reflect_tidy_corr(tidy_corr_all)
    tidy_violin_all['cellpop'] = 'all'

    tidy_violin_resp = reflect_tidy_corr(tidy_corr_resp)
    tidy_violin_resp['cellpop'] = 'responders'

    tidy_violin = tidy_violin_all.append(tidy_violin_resp)

    tidy_to_plot = tidy_violin.loc[tidy_violin['corrtype'] == corrtype, :]
    g = sns.catplot(data=tidy_to_plot, x="stimname_col", y="value", hue="cellpop",
                    col="stimname_row", col_wrap=3,
                    kind="violin", height=4.5, aspect=0.75,
                    split=True, inner='stick',
                    palette="pastel",
                    tight_layout=True, legend_out=True)
    return g


def dist_to_tidy_tril_corr(dist):
    """Converts dict of correlation matrices into """
    df_list = []
    for k in dist.keys():
        df = tidy.tidy_tril_df_corrmat(dist[k])
        df['corrtype'] = k
        df_list.append(df)
    tidy_corr = pd.concat(df_list, axis=0)  # lower triangular correlation info
    return tidy_corr


def dist_to_tidy_triu_corr(dist):
    df_list = []
    for k in dist.keys():
        df = tidy.tidy_triu_df_corrmat(dist[k])
        df['corrtype'] = k
        df_list.append(df)
    tidy_corr = pd.concat(df_list, axis=0)  # lower triangular correlation info
    return tidy_corr


def extract(f):
    """
    Load traces from suite2p/**/step00_process_suite2p (folder in the same directory as stat.npy)

    :param f: Path to stat.npy, uses this to parse fly information and load other dataset info
    :return:
    """
    nt_meta_data = postprocess_metadata(load_metadata(f))
    ainfo = dirutils.parse_dataset_info(f)

    title_str = get_title_str(ainfo, nt_meta_data)

    folder = f.parent
    C_dec = np.load(folder.joinpath('step00_process_suite2p', 'C_dec.npy'))
    return nt_meta_data, title_str, C_dec


def main(f):
    # %%
    ainfo = dirutils.parse_dataset_info(f)
    pp.pprint(ainfo)

    # loading metadata
    if ainfo.is_block:
        meta_data_loader = analysis.BlockMetaDataLoader.init_from_statfile(f)
    else:
        meta_data_loader = analysis.MovieMetaDataLoader.init_from_statfile(f)
    nt_meta_data = postprocess_metadata(meta_data_loader.load())

    # loading suite2p output data
    if f.parent.name == 'combined':
        suite2p_loader = analysis.Suite2pCombinedLoader(statfile=f)
    elif 'plane' in f.parent.name:
        suite2p_loader = analysis.Suite2pPlane0Loader(statfile=f)
    nt_s2p_data = suite2p_loader.load()

    # %% save nt_meta_data to .mat file
    nt_meta_data_dict = dict(nt_meta_data._asdict())
    nt_meta_data_dict['df_frametimeinfo'] = nt_meta_data_dict['df_frametimeinfo'].to_records(index=False)
    nt_meta_data_dict['df_stimulus'] = nt_meta_data_dict['df_stimulus'].to_records(index=False)
    nt_meta_data_dict['meta_yaml'] = json.dumps(nt_meta_data.meta_yaml)
    savemat(f.with_name('nt_meta_data.mat'), nt_meta_data_dict)

    # compute frame rates
    print('computing frame rates')
    fs = trace.compute_fs(nt_meta_data.timestamps)
    print('done')

    # extract loaded files
    df_stimulus = copy.deepcopy(nt_meta_data.df_stimulus)
    df_frametimeinfo = copy.deepcopy(nt_meta_data.df_frametimeinfo)
    timestamps = copy.deepcopy(nt_meta_data.timestamps)
    ti = copy.deepcopy(nt_meta_data.sync_meta)
    iscell = nt_s2p_data.iscell[:, 0].astype(np.bool_)
    cellprob = nt_s2p_data.iscell[:, 1]

    use_saved_cells = True
    if use_saved_cells:
        good_cells_bool = iscell
    else:
        good_cells_bool = (cellprob > 0.5).squeeze()
    good_cells = np.argwhere(good_cells_bool).squeeze()

    # save trial tensors of raw fluorescence traces
    tensor_F, trial_ts = represp.align_signals(nt_s2p_data.F, timestamps, df_stimulus['stim_ict'])
    tensor_Fneu, trial_ts = represp.align_signals(nt_s2p_data.Fneu, timestamps, df_stimulus['stim_ict'])

    np.savez(f.with_name("TENSOR_RAW.npz"), tensor_F=tensor_F, tensor_Fneu=tensor_Fneu, trial_ts=trial_ts)
    savemat(f.with_name("TENSOR_RAW.mat"), dict(tensor_F=tensor_F, tensor_Fneu=tensor_Fneu, trial_ts=trial_ts))

    # compute block concentration
    blk_conc = df_stimulus.blk_conc.unique()[0]

    # rearing condition
    rear_type = nt_meta_data.session_info['rearing']
    # %% make title for plots
    title_str = f"""{ainfo.flydate} fly {ainfo.flynum}, {ainfo.mov}"""
    if ainfo.is_block:
        title_str = f"""{title_str} blk={ainfo.blk}\n"""

    title_str = title_str + f"""\nreared: {rear_type}    |   conc={blk_conc}"""
    print(title_str)

    # %% calculate zscored/dfof traces, don't drop cells yet
    print('processing fluorescent traces')
    F = nt_s2p_data.F - 0.7 * nt_s2p_data.Fneu
    F_blocked = trace.split_traces_to_blocks(F, nt_meta_data.df_frametimeinfo.scope_pulse_idx)

    # BASELINE
    medfilt = partial(filters.median_filter, size=(1, 4))
    gauss1 = partial(gaussian_filter1d, sigma=1)
    gauss1s = partial(gaussian_filter1d, sigma=int(fs))
    gauss3s = partial(gaussian_filter1d, sigma=int(3 * fs))
    gauss10s = partial(filters.gaussian_filter1d, sigma=10 * int(fs))

    # zero w/ 5th percentile, then zscore
    low = [np.percentile(item, 5, axis=1, keepdims=True) for item in F_blocked]
    Fz = [zscore(a - b, axis=1) for a, b in zip(F_blocked, low)]

    # lowpass filter, smooth, rolling percentile, smooth, smooth
    Ffilt = [lowpass_filter(item, Fs=fs, cutoff=0.5) for item in Fz]
    Fsmooth = [gauss1s(item) for item in Ffilt]
    Flow = [filters.percentile_filter(item, percentile=5, size=(1, int(90 * fs))) for item in Fsmooth]
    Flow = [gauss10s(gauss10s(item)) for item in Flow]

    # calculate rolling percentile baseline
    dfof = [a - b for a, b in zip(Ffilt, Flow)]
    zdfof = [zscore(item) for item in dfof]
    dfof = np.hstack(tuple(dfof))
    zdfof = np.hstack(tuple(zdfof))

    save_zdfof = True
    if save_zdfof:
        np.savez(f.with_name('ZDFOF.npz'), dfof=dfof, zdfof=zdfof)
        savemat(f.with_name('ZDFOF.mat'), dict(dfof=dfof, zdfof=zdfof))
    # %% TRIAL X CELL X TIME TENSORS

    zdfof_tensor, trial_ts = represp.align_signals(zdfof, timestamps, df_stimulus['stim_ict'])
    dfof_tensor, trial_ts = represp.align_signals(dfof, timestamps, df_stimulus['stim_ict'])
    ustim = df_stimulus.stimname.unique().categories

    # # calculate means across odors
    # mean_odors = []
    # for iodor, odor in enumerate(ustim):
    # 	mean_odors.append(np.mean(zdfof_tensor[:, df_stimulus.stimname == odor, :], axis=1))
    #
    # mean_odors_zdfof = np.stack(mean_odors)
    # mean_zdfof_hcat = np.hstack(mean_odors)

    TENSOR_TRIALS = dict(zdfof_tensor=zdfof_tensor, dfof_tensor=dfof_tensor,
                         trial_ts=trial_ts, stimname=df_stimulus.stimname.to_list(),
                         stim_idx=df_stimulus.stim_idx.to_list())
    np.savez(f.with_name("TENSOR_TRIALS.npz"), **TENSOR_TRIALS)
    savemat(f.with_name("TENSOR_TRIALS.mat"), TENSOR_TRIALS)

    # %% STACKED PLOT
    cid0 = 0
    cid_idx = np.arange(20) + cid0
    cids = list(good_cells[cid_idx])
    fig, axarr = vis.blocked_stackplot(dfof, timestamps, cids=cids,
                                       df_stimulus=df_stimulus,
                                       df_frametimeinfo=df_frametimeinfo,
                                       peak_start=1, peak_len=5)
    plt.show()
    fig.savefig(f.with_name('blocked_stackedplot.pdf'))
    # fig.suptitle(f"""{title_str}\n{json.dumps(dfof_s2p_ops)} """, fontsize=10)
    # fig.suptitle(title_str + f"""\n\nzdfof settings: {json.dumps(dfof_s2p_ops)}""")
    # plt.show()

    # %% compute response vectors

    respvec_ops = dict(peak_start=1, peak_len=5, baseline_len=-10)
    stim_ict = ti['stim_ict'][df_stimulus['stim_idx'].to_numpy() - 1]
    df_stimulus['stim_ict'] = stim_ict

    peak_mean, baseline_mean, baseline_std = respvec.compute_interval_response_vectors(zdfof,
                                                                                       timestamps,
                                                                                       df_stimulus['stim_ict'],
                                                                                       **respvec_ops)
    peak_amp = peak_mean - baseline_mean
    peak_nstd = peak_amp / baseline_std

    save_respvec = True
    if save_respvec:
        RESPVEC = dict(peak_amp=peak_amp, peak_nstd=peak_nstd, baseline_std=baseline_std, baseline_mean=baseline_mean,
                       respvec_ops=respvec_ops)
        np.savez(f.with_name('RESPVEC.npz'), **RESPVEC)
        savemat(f.with_name('RESPVEC.mat'), RESPVEC)
    # np.save(f.with_name('respvec_peak_amp.npy'), peak_amp)
    # np.save(f.with_name('respvec_peak_nstd.npy'), peak_nstd)
    # np.save(f.with_name('respvec_baseline_std.npy'), baseline_std)
    # np.save(f.with_name('respvec_baseline_mean.npy'), baseline_mean)
    # np.save(f.with_name('respvec_ops.npy'), respvec_ops)
    # %% compute responders
    bin_response = respvec.get_binresponse(peak_amp, baseline_std, n_std=3, thr=0.5)
    bin_idx = np.any(bin_response, axis=1)

    # response fraction plot

    response_frac = np.sum(bin_response[good_cells, :], axis=0) / good_cells.size
    fig, ax = plt.subplots(figsize=(4, 4), tight_layout=True)
    ax = sns.barplot(x=df_stimulus['stimname'][1:], y=response_frac[1:], palette=sns.color_palette("hls", 8))
    ax.set_title(title_str + f"""\nResponse fraction""", fontsize=10)
    ax.set_ylabel('fraction responders')
    ax.set_ylim([0, 0.25])
    plt.show()

    save_sparsity_plot = True
    if save_sparsity_plot:
        fig.savefig(f.with_name('response_strength.pdf'))
        fig.savefig(f.with_name('response_strength.png'))

    # %%  plot correlation matrices
    n_cells = good_cells.size
    n_responder_cells = np.sum(bin_idx & good_cells_bool)

    # compute correlation matrices
    dist_all, df_peakresp = respvec.compute_corrmats(peak_amp[good_cells, :], df_stimulus)
    dist_resp, df_peakresp = respvec.compute_corrmats(peak_amp[bin_idx & good_cells_bool, :], df_stimulus)
    dist_all_mean = mean_corrmats(dist_all)
    dist_resp_mean = mean_corrmats(dist_resp)
    np.savez(f.with_name("DIST.npz"), dist_all=dist_all, dist_all_mean=dist_all_mean,
             dist_resp=dist_resp, dist_resp_mean=dist_resp_mean)

    with PdfPages(f.with_name('corrmat_multipage.pdf')) as pdf:
        fig_corrmat_all, axarr, haxarr = vis.corrmats(dist_all)
        fig_corrmat_all.suptitle(f"""{title_str}\nall cells (n_cells={n_cells})""")
        pdf.savefig(fig_corrmat_all)
        plt.show()

        # fig_corrmat_allmean
        # ----------------------------------------------------------------------
        fig_corrmat_allmean, axarr, haxarr = vis.mean_corrmats(dist_all_mean, annot=True)
        fig_corrmat_allmean.suptitle(f"""{title_str}\nall cells (n_cells={n_cells})""")
        pdf.savefig(fig_corrmat_allmean)
        plt.show()

        # fig_corrmat_resp: plot correlation matrices, with only responders included
        # ----------------------------------------------------------------------
        fig_corrmat_resp, axarr, haxarr = vis.corrmats(dist_resp)
        fig_corrmat_resp.suptitle(f"""{title_str}\nresponders only ({n_responder_cells}/{n_cells})""")
        pdf.savefig(fig_corrmat_resp)
        plt.show()

        # fig_corrmat_respmean
        # ----------------------------------------------------------------------
        fig_corrmat_respmean, axarr, haxarr = vis.mean_corrmats(dist_resp_mean, annot=True)
        fig_corrmat_respmean.suptitle(f"""{title_str}\nresp. only ({n_responder_cells}/{n_cells})""")
        pdf.savefig(fig_corrmat_respmean)
        plt.show()

    print('multipage corrmat plot pdf saved.')
    # ----------------------------------------------------------------------

    # %%
    save_pdf_plots = False
    if save_pdf_plots:
        print('saving corrmat plots as pdf')
        fig_corrmat_all.savefig(f.with_name('corrmat_all.pdf'))
        fig_corrmat_allmean.savefig(f.with_name('corrmat_all_mean.pdf'))

        fig_corrmat_resp.savefig(f.with_name('corrmat_resp.pdf'))
        fig_corrmat_respmean.savefig(f.with_name('corrmat_resp_mean.pdf'))
        print('corrmat plots saved as pdf')

    save_png_plots = False
    if save_png_plots:
        print('saving corrmat plots as png')
        fig_corrmat_all.savefig(f.with_name('corrmat_all.png'))
        fig_corrmat_allmean.savefig(f.with_name('corrmat_all_mean.png'))

        fig_corrmat_resp.savefig(f.with_name('corrmat_resp.png'))
        fig_corrmat_respmean.savefig(f.with_name('corrmat_resp_mean.png'))
        print('corrmat plots saved as png')

    # %%
    def add_info_to_tidy_corr(df_):
        """ Adds relevant dataset info (ainfo fields) to tidy correlation info plots"""
        df_.reset_index()
        df_.reset_index(drop=True, inplace=True)
        df_['flydate'] = ainfo.flydate
        df_['flynum'] = ainfo.flynum
        df_['mov'] = ainfo.mov
        df_['rear_type'] = rear_type
        if ainfo.is_block:
            df_['chunk'] = ainfo.blk
        else:
            df_['chunk'] = np.nan
        df_['conc'] = blk_conc
        return df_

    tidy_corr_all_triu = dist_to_tidy_triu_corr(dist_all)
    tidy_corr_all_triu = add_info_to_tidy_corr(tidy_corr_all_triu)
    tidy_corr_all_triu['cellpop'] = 'all'

    tidy_corr_all_tril = dist_to_tidy_tril_corr(dist_all)
    tidy_corr_all_tril = add_info_to_tidy_corr(tidy_corr_all_tril)
    tidy_corr_all_tril['cellpop'] = 'all'

    tidy_corr_resp_triu = dist_to_tidy_triu_corr(dist_resp)
    tidy_corr_resp_triu = add_info_to_tidy_corr(tidy_corr_resp_triu)
    tidy_corr_resp_triu['cellpop'] = 'responders'

    tidy_corr_resp_tril = dist_to_tidy_triu_corr(dist_resp)
    tidy_corr_resp_tril = add_info_to_tidy_corr(tidy_corr_resp_tril)
    tidy_corr_resp_tril['cellpop'] = 'responders'

    tidy_corr_triu = tidy_corr_all_triu.append(tidy_corr_resp_triu)
    tidy_corr_tril = tidy_corr_all_tril.append(tidy_corr_resp_tril)

    save_tidy_corr = True
    if save_tidy_corr:
        feather.write_feather(tidy_corr_all_triu, f.with_name("tidy_corr_all_triu.feather"), compression='uncompressed')
        print('saved: tidy_corr_all_triu.feather')
        feather.write_feather(tidy_corr_all_tril, f.with_name("tidy_corr_all_tril.feather"), compression='uncompressed')
        print('saved: tidy_corr_all_tril.feather')
        feather.write_feather(tidy_corr_resp_triu, f.with_name('tidy_corr_resp_triu.feather'),
                              compression='uncompressed')
        print('saved: tidy_corr_resp_triu.feather')
        feather.write_feather(tidy_corr_resp_tril, f.with_name('tidy_corr_resp_tril.feather'),
                              compression='uncompressed')
        print('saved: tidy_corr_resp_tril.feather')
        feather.write_feather(tidy_corr_tril, f.with_name('tidy_corr_tril.feather'),
                              compression='uncompressed')
        print('saved: tidy_corr_tril.feather')
        feather.write_feather(tidy_corr_triu, f.with_name('tidy_corr_triu.feather'),
                              compression='uncompressed')
        print('saved: tidy_corr_triu.feather')
    # %%
    with PdfPages(f.with_name('split_violin_corrscatter_plots.pdf')) as pdf:
        for corrtype in ['cosine', 'pearson', 'kendall', 'spearman']:
            tidy_to_plot = reflect_tidy_corr(tidy_corr_triu.loc[tidy_corr_triu['corrtype'] == corrtype, :])
            g = sns.catplot(data=tidy_to_plot, x="stimname_col", y="value", hue="cellpop",
                            col="stimname_row", col_wrap=3,
                            kind="violin", height=4.5, aspect=0.6,
                            split=True, inner='stick',
                            palette="pastel",
                            tight_layout=True, legend_out=True)
            plt.subplots_adjust(top=0.9, left=0.05, bottom=0.1)  # adjust spacing
            g.fig.suptitle(f"""{title_str}\ndistance metric={corrtype} corr.""", fontsize=10)
            g.set_axis_labels("", corrtype)
            g.set(ylim=(-0.4, 1.0))
            g.map(plt.axhline, y=0, linewidth=0.5, color='k')
            pdf.savefig(g.fig)
            plt.show()
    print(f"saved {f.with_name('split_violin_corrscatter_plots.pdf')}")


# %%
if __name__ == '__main__':
    stat_files = [Path(r".\data\processed_data\2021-07-27\1\movie_003\suite2p\combined\stat.npy")]
    for f in stat_files:
        main(f)

# # odor reared
# odor_files = [Path(r".\data\processed_data\2021-06-26\1\movie_005\suite2p\plane0\stat.npy")]
#
# # pfo reared
# pfo_files = [Path(r".\data\processed_data\2021-06-26\2\movie_001\suite2p\plane0\stat.npy")]
#
# stat_files = [*pfo_files, *odor_files]

# %%


#
#
# # %%
#
# def main(f):
#     """
#     ** Note: rear_type calculated here **
#
#     Parse the filepath to ** embedding.npy ** to extract info about  the  dataset, including rearing condition.
#         - `flydate`
#         - `flynum`
#         - `mov`
#         - `blk`
#
#     Load analysis outputs  from the suite2p folder.
#         - `Fcell` -from F.npy
#         - `Fneu` - from Fneu.npy
#         - `iscell` - from iscell.npy (first column = 0 or 1, second column = probability of being a cell)
#         - `stat` - from stat.npy
#         - `embedding` - from embedding.npy
#
#     Depending on if the dataset 's movie was split into chunks (sub-movies) for analysis,
#     load `df_frametimeinfo`, `df_stimulus`, and `timestamps` from the chunk directory (`blk_dir` )
#     or the movie directory (`blk_dir`)
#
#     Parameters
#     ----------
#     f : pathlib.Path
#     """
#     # get movie directory
#
#     ainfo = dirutils.parse_filepath(f)  # analysis info
#
#     mov_dir = parse_parent_movie(f)
#
#     # get experiment info
#     flydate = parse_parent_date(f).name
#     flynum = parse_parent_fly(f).name
#     mov = parse_parent_movie(f).name
#
#     # load session info, determine rearing type, load stimulus timing info from sync_meta.json
#     session_info = pd.read_excel(mov_dir.joinpath('metadata.xlsx'), sheet_name='sessioninfo').squeeze()
#     rear_type = session_info['rearing']
#     ti = thorsync.load_sync_meta(mov_dir.joinpath("sync_meta.json"))
#
#     if ainfo.is_block:
#         blk_dir = parse_parent_block(f)
#         blk = parse_parent_block(f).name
#
#         # load df_stimulus, df_frametimeinfo, and timestamps from blk_dir
#         df_stimulus = pd.read_csv(blk_dir.joinpath('metadata.csv'))
#         df_frametimeinfo = pd.read_csv(blk_dir.joinpath('frametimeinfo.csv'))
#         timestamps = np.load(blk_dir.joinpath('timestamps.npy'))
#         title_str = f"""{flydate} fly {flynum}, {mov} block {blk}, \n{rear_type},"""
#
#     else:
#         # load df_stimulus, df_frametimeinfo, and timestamps from mov_dir
#         df_stimulus = pd.read_excel(mov_dir.joinpath('metadata.xlsx'), sheet_name='stimulus')
#         df_frametimeinfo = pd.read_csv(mov_dir.joinpath('frametimeinfo.csv'))
#         timestamps = np.load(mov_dir.joinpath('timestamps.npy'))
#         title_str = f"""{flydate} fly {flynum}, {mov}\n{rear_type},"""
#
#     # load suite2p output files
#     Fcell = np.load(f.with_name('F.npy'))
#     Fneu = np.load(f.with_name("Fneu.npy"))
#     iscell = np.load(f.with_name("iscell.npy"))
#     stat = np.load(f.with_name("stat.npy"), allow_pickle=True)
#     embedding = np.load(f.with_name('embedding.npy'), allow_pickle=True).item()
#
#     # subtract background neuropil activity
#     F = Fcell - 0.7 * Fneu
#     cellprob = iscell[:, 1]
#     iscell = iscell[:, 0].astype(np.bool_)
#
#     F = F[iscell, :]
#     n_cells = F.shape[0]
#
#     # compute sampling rate from timestamps
#     fs = trace.compute_fs(timestamps)
#
#     n_trials = df_stimulus.shape[0]
#
#     # compute block/chunk concentration, from df_stimulus [ideally there's only 1 stimulus type per movie]
#     blk_conc = min(df_stimulus.loc[:, ['conc_a', 'conc_b']].mode(axis=0).min(axis=0).tolist())
#
#
#
#     # normalize/preprocess neuropil-corrected fluorescence traces
#     sp = copy.deepcopy(F)
#     spn = rasterplot.norm_traces(sp)
#     spp = rasterplot.prep_to_plot(spn)
#
#     print(f"\nfile: {f}")
#     print("----------------------------------------------------------------")
#     print(f"\n{title_str}")
#
#     # ------------------------------------------------------------------------------------------------------------------
#     # Correct df_frametimeinfo
#     # ------------------------------------------------------------------------------------------------------------------
#     # If df_frametimeinfo contains timestamps for multiple planes, drop everything but plane=0 timestamps,
#     # and check that the dimensions match timestamps
#
#     if 'plane' in df_frametimeinfo.columns:
#         df_frametimeinfo = df_frametimeinfo.loc[df_frametimeinfo['plane'] == 0]
#         df_frametimeinfo.reset_index(drop=True, inplace=True)  #
#
#     print(f"\nsize of timestamps array: {timestamps.size}")
#     print(f"\nsize of df_frametimeinfo dataframe: {df_frametimeinfo.shape}")
#
#     if timestamps.size != df_frametimeinfo.shape[0]:
#         raise ValueError("Dimensions of timestamps and df_frametimeinfo are not compatible! ")
#
#     # ------------------------------------------------------------------------------------------------------------------
#     # Add info to df_stimulus
#     # ------------------------------------------------------------------------------------------------------------------
#     # calculate timing of stimuli, using ti['stim_ict']
#     df_stimulus['stim_ict'] = ti['stim_ict'][df_stimulus['stim_idx'].to_numpy() - 1]
#
#     df_stimulus['scope_pulse_idx'] = np.interp(df_stimulus['stim_ict'],
#                                                df_frametimeinfo['frame_times'],
#                                                df_frametimeinfo['scope_pulse_idx']).astype('int')
#     df_stimulus['stimname'] = df_stimulus['stimname'].astype(ocat)
#
#     # print(df_stimulus)
#     # ------------------------------------------------------------------------------------------------------------------
#     # split traces and lowpass filter
#     # ------------------------------------------------------------------------------------------------------------------
#     blocked_traces = trace.split_traces_to_blocks(F, df_frametimeinfo.scope_pulse_idx)
#
#     filt_ops = dict(fs=fs, cutoff=0.5, trans_width=0.05, numtaps=230)
#     b = signal.remez(filt_ops['numtaps'],
#                      [0, filt_ops['cutoff'],
#                       filt_ops['cutoff'] + filt_ops['trans_width'],
#                       0.5 * filt_ops['fs']], [1, 0], Hz=filt_ops['fs'])
#     filt_blocked_traces = [signal.filtfilt(b, 1, y) for y in blocked_traces]
#
#     # ------------------------------------------------------------------------------------------------------------------
#     # compute dfof, in one of two ways
#     # ------------------------------------------------------------------------------------------------------------------
#     dfof_method = 'trace'
#
#     # suite2p preprocessing
#     if dfof_method == 'suite2p':
#         blocked_dfof, blocked_df, blocked_baseline = zip(*[trace.suite2p_extract_dfof(
#             F=item,
#             baseline='constant_prctile',  # 'constant_prctile', or 'maximin'
#             win_baseline=90,
#             sig_baseline=5,
#             fs=fs,
#             prctile_baseline=10) for item in filt_blocked_traces])
#
#     # rolling percentile window
#     if dfof_method == 'trace':
#         blocked_dfof, blocked_df, blocked_baseline = zip(*[trace.extract_dfof(item,
#                                                                               fs=fs,
#                                                                               win=60,
#                                                                               percentile=10)
#                                                            for item in blocked_traces])
#
#     blocked_zdfof = [zscore(item, axis=1) for item in filt_blocked_traces]
#
#     # combine lists of traces
#     dfof = np.hstack(tuple(blocked_dfof)).astype('float32')
#     zdfof = np.hstack(tuple(blocked_zdfof)).astype('float32')
#     np.save(f.with_name('dfof.npy'), dfof)
#     np.save(f.with_name('zdfof.npy'), zdfof)
#
#     # dfof_tensor = represp.align_signals(dfof, timestamps, df_stimulus['stim_ict'], response_len=30, baseline_len=15)
#     # zdfof_tensor = represp.align_signals(dfof, timestamps, df_stimulus['stim_ict'], response_len=30, baseline_len=15)
#     # np.save(f.with_name('dfof_tensor.npy'), dfof)
#     # np.save(f.with_name('zdfof_tensor.npy'), dfof)
#     # ------------------------------------------------------------------------------------------------------------------
#     # compute peak response amplitudes
#     # ------------------------------------------------------------------------------------------------------------------
#
#     # compute response vectors
#     # peak_mean, baseline_mean, baseline_std = respvec.compute_interval_response_vectors(zdfof,
#     #                                                                                    timestamps,
#     #                                                                                    ti['stim_ict'],
#     #                                                                                    peak_start=1,
#     #                                                                                    peak_len=5,
#     #                                                                                    baseline_len=-10)
#     # %%
#
#     respvec_ops = dict(peak_start=1, peak_len=5, baseline_len=-10)
#     stim_ict = ti['stim_ict'][df_stimulus['stim_idx'].to_numpy() - 1]
#     df_stimulus['stim_ict'] = stim_ict
#
#     peak_mean, baseline_mean, baseline_std = respvec.compute_interval_response_vectors(zdfof,
#                                                                                        timestamps,
#                                                                                        df_stimulus['stim_ict'],
#                                                                                        **respvec_ops)
#     peak_amp = peak_mean - baseline_mean
#     peak_nstd = peak_amp / baseline_std
#
#     # compute responders
#     bin_response = respvec.get_binresponse(peak_amp, baseline_std, n_std=3, thr=0.5)
#     bin_idx = np.any(bin_response, axis=1)
#
#     # compute correlation matrices
#     dist_all, df_peakresp = respvec.compute_corrmats(peak_amp, df_stimulus)
#     dist_resp, df_peakresp = respvec.compute_corrmats(peak_amp[bin_idx, :], df_stimulus)
#     dist_all_mean = mean_corrmats(dist_all)
#     dist_resp_mean = mean_corrmats(dist_resp)
#
#     # fig corrmat_all: plot correlation matrices, with all cells included
#     fig_corrmat_all, axarr, haxarr = vis.corrmats(dist_all)
#     fig_corrmat_all.suptitle(f"""{title_str}\nall cells""")
#     plt.show()
#
#     # fig_corrmat_allmean
#     fig_corrmat_allmean, axarr, haxarr = vis.mean_corrmats(dist_all_mean, annot=True)
#     fig_corrmat_allmean.suptitle(f"""{title_str}\nall cells (mean)""")
#     plt.show()
#
#     # fig_corrmat_resp: plot correlation matrices, with only responders included
#     fig_corrmat_resp, axarr, haxarr = vis.corrmats(dist_resp)
#     fig_corrmat_resp.suptitle(f"""{title_str}\nresponders only""")
#     plt.show()
#
#     # fig_corrmat_respmean
#     fig_corrmat_respmean, axarr, haxarr = vis.mean_corrmats(dist_resp_mean, annot=True)
#     fig_corrmat_respmean.suptitle(f"""{title_str}\nresp. only (mean)""")
#     plt.show()
#
#     #%% save corrmat plots
#
#     save_pdf_plots = True
#     if save_pdf_plots:
#         print('saving corrmat plots as pdf')
#         fig_corrmat_all.savefig(f.with_name('corrmat_all.pdf'))
#         fig_corrmat_allmean.savefig(f.with_name('corrmat_all_mean.pdf'))
#
#         fig_corrmat_resp.savefig(f.with_name('corrmat_resp.pdf'))
#         fig_corrmat_respmean.savefig(f.with_name('corrmat_resp_mean.pdf'))
#         print('corrmat plots saved as pdf')
#
#     save_png_plots = True
#     if save_png_plots:
#         print('saving corrmat plots as png')
#         fig_corrmat_all.savefig(f.with_name('corrmat_all.png'))
#         fig_corrmat_allmean.savefig(f.with_name('corrmat_all_mean.png'))
#
#         fig_corrmat_resp.savefig(f.with_name('corrmat_resp.png'))
#         fig_corrmat_respmean.savefig(f.with_name('corrmat_resp_mean.png'))
#         print('corrmat plots saved as png')
#
#
#     # %% save dist files
#
#     save_dist=True
#     if save_dist:
#         np.save(f.with_name('dist_all.npy'), dist_all, allow_pickle=True)
#         np.save(f.with_name('dist_resp.npy'), dist_resp, allow_pickle=True)
#         np.save(f.with_name('dist_all_mean.npy'), dist_all_mean, allow_pickle=True)
#         np.save(f.with_name('dist_resp_mean.npy'), dist_resp_mean, allow_pickle=True)
#
#     # %% accumulate and save correlations from correlation matrices
#
#
#     #%%
#
#
# # %%
# if __name__ == '__main__':
#
#     DATA_DIR = Path(r'./data/processed_data')
#
#     file_list = sorted(list(DATA_DIR.rglob(r"2021-06-26/1/**/suite2p/plane0/embedding.npy")))
#
#     # print files
#     for i, file in enumerate(file_list):
#         print(f"{i} ----- {file}")
#
#     # ask user to select which files to run
#     file_idx = int(input("Enter file #: (enter -1 to run all)"))
#     if file_idx == -1:
#         print("RUNNING ALL: ")
#         print("------------")
#     else:
#         file = file_list[file_idx]
#
#     # loop over file list, and run correlation analysis
#     for file in file_list:
#         dist_all, dist_resp, dist_all_mean, dist_resp_mean, title_str = main(file)
#
#         with open(file.with_name('title_str.txt'), "w") as text_file:
#             text_file.write(title_str)
#
#         # make violin plots
#         for corrtype in ['cosine', 'pearson', 'spearman', 'kendall']:
#             g = faceted_corrscatter(file, corrtype)
#
#             plt.subplots_adjust(top=0.9, left=0.05, bottom=0.1)  # adjust spacing
#             g.fig.suptitle(f"""{title_str}\ndistance metric={corrtype} corr.""", fontsize=10)
#             g.set_axis_labels("", corrtype)
#             plt.show()
#
#             # save plots to pdf
#             g.fig.savefig(file.with_name(f"violin_corrscatter__{corrtype}.png"))
#             g.fig.savefig(file.with_name(f"violin_corrscatter__{corrtype}.pdf"))
#
#         # dist_all, dist_resp, dist_all_mean, dist_resp_mean, title_str = main(f)
