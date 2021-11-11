import copy
import importlib as imp
import json
from collections import namedtuple
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.io as pio
from scipy import signal
from scipy.stats import zscore

pio.renderers.default = "browser"

from suite2p.extraction import dcnv

from rearing.ca import trace, vis, represp
from ryim import dirutils, analysis
import rasterplot
from rastermap import Rastermap

import pprint as pp

# import nitime
# import nitime.timeseries as ts
# import nitime.analysis as nta
# import nitime.viz as viz

# %%
ocat = pd.CategoricalDtype(categories=['pfo', '1-5ol', 'VA', '1-6ol', 'EP', '1-6ol+EP'], ordered=True)

MetaData = namedtuple("MetaData", ["df_stimulus", "df_frametimeinfo", "timestamps",
                                   "session_info", "xml_meta", "sync_meta", "meta_yaml"])


# %%
def postprocess_metadata(_nt_meta_data, ocat=ocat):
    """
    fix some of the metadata stuff
    :param _nt_meta_data:
    :param ocat:
    :return:
    """

    _ti = _nt_meta_data.sync_meta
    _df_stimulus = _nt_meta_data.df_stimulus
    _df_frametimeinfo = _nt_meta_data.df_frametimeinfo

    _df_stimulus['blk_conc'] = min(_df_stimulus.loc[:, ['conc_a', 'conc_b']].mode(axis=0).min(axis=0).tolist())
    _df_stimulus['stim_ict'] = _ti['stim_ict'][_df_stimulus['stim_idx'].to_numpy() - 1]
    _df_stimulus['scope_pulse_idx'] = np.interp(_df_stimulus['stim_ict'],
                                                _df_frametimeinfo['frame_times'],
                                                _df_frametimeinfo['scope_pulse_idx']).astype('int')

    _df_stimulus['stimname'] = _df_stimulus['stimname'].astype(ocat)

    return MetaData(_df_stimulus,
                    _df_frametimeinfo,
                    _nt_meta_data.timestamps,
                    _nt_meta_data.session_info,
                    _nt_meta_data.xml_meta,
                    _nt_meta_data.sync_meta,
                    _nt_meta_data.meta_yaml
                    )


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


def compute_blocked_suite2p_Fc(_list_of_cells_x_time, Fs, **kwargs):
    """
    Use suite2p's df/f function.

    Default parameters are dict(baseline='maximin', win_baseline=60, sig_baseline=5, prctile_baseline=10)

    :param _list_of_cells_x_time:
    :param Fs:
    :param kwargs:
    :return:
        list blk_dfof
        list blk_df
        list blk
    """

    dfof_ops = dict(baseline='maximin', win_baseline=60, sig_baseline=5, prctile_baseline=10)
    dfof_ops.update(kwargs)
    print(dfof_ops)

    nblocks = len(_list_of_cells_x_time)
    blk_Fc = [None] * nblocks
    blk_df = [None] * nblocks
    blk_baseline = [None] * nblocks



    for i, _cells_x_time in enumerate(_list_of_cells_x_time):
        _dfof, _df, _baseline = F_corrected = dcnv.preprocess(F=F, fs=Fs, **dfof_ops)
        blk_dfof[i] = _dfof
        blk_df[i] = _df
        blk_baseline[i] = _baseline
    return blk_dfof, blk_df, blk_baseline


def compute_blocked_prctile_dfof(_list_of_cells_x_time, Fs, **kwargs):
    """
    Compute df/f with rolling percentile method.
    Default parameters are:
        win: 60 (sec)
        percentile: 10

    :param _list_of_cells_x_time:
    :param Fs:
    :param kwargs:
    :return:
    """
    dfof_ops = dict(win=60, percentile=10)
    dfof_ops.update(kwargs)

    n_blocks = len(_list_of_cells_x_time)
    blk_dfof = [None] * n_blocks
    blk_df = [None] * n_blocks
    blk_baseline = [None] * n_blocks

    for i, _cells_x_time in enumerate(_list_of_cells_x_time):
        _dfof, _df, _baseline = trace.extract_dfof(_cells_x_time, Fs, **dfof_ops)
        blk_dfof[i] = _dfof
        blk_df[i] = _df
        blk_baseline[i] = _baseline
    return blk_dfof, blk_df, blk_baseline


# %%
# set some default parameters

# default suite2p dfof params
dfof_s2p_ops = dict(baseline='maximin', win_baseline=60, sig_baseline=5, prctile_baseline=20)

# default dfof params
dfof_prctile_ops = dict(win=60, percentile=10)

# default rasterplot params
rasterplot_ops = {'n_components': 1, 'n_X': 100, 'alpha': 1., 'K': 1.,
                  'nPC': 200, 'constraints': 2, 'annealing': True, 'init': 'pca',
                  }
# %%

#stat_file_1 = Path(r".\data\processed_data\2021-08-30\4\movie_001\suite2p\combined\stat.npy")
#statfiles = [stat_file_1]
#stat_file_1 = Path("./data/processed_data/2021-07-27/1/movie_003/suite2p/combined/stat.npy")

#f = Path(r".\data\processed_data\2021-08-30\2\movie_001\1\suite2p\combined\stat.npy")
#statfiles = [stat_file_1]

#%%
def main(f):
    #%%
    # ainfo = AnalysisSet(flydate, flynum, mov, is_blk, blk)
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

    # compute frame rates
    print('computing frame rates')
    fs = trace.compute_fs(nt_meta_data.timestamps)
    print('done')

    # extract loaded files
    df_stimulus = copy.deepcopy(nt_meta_data.df_stimulus)
    df_frametimeinfo = copy.deepcopy(nt_meta_data.df_frametimeinfo)
    timestamps = copy.deepcopy(nt_meta_data.timestamps)
    iscell = nt_s2p_data.iscell[:, 0].astype(np.bool_)
    ti = copy.deepcopy(nt_meta_data.sync_meta)

    # compute block concentration
    blk_conc = df_stimulus.blk_conc.unique()[0]

    # rearing condition
    rear_type = nt_meta_data.session_info.rearing
    # -----------------------------------------------
    # process fluorescence traces
    # -----------------------------------------------
#%%	# %% make title for plots
    title_str = f"""{ainfo.flydate} fly {ainfo.flynum}, {ainfo.mov}"""

    if ainfo.is_block:
        title_str = f"""{title_str} blk={ainfo.blk}\n"""

    title_str = title_str + f"""\nreared: {rear_type}    |   conc={blk_conc}"""
    print(title_str)

    # %% calculate zscored/dfof traces
    print('processing fluorescent traces')

    F = nt_s2p_data.F - 0.7 * nt_s2p_data.Fneu

    F_blocked = trace.split_traces_to_blocks(F, nt_meta_data.df_frametimeinfo.scope_pulse_idx)
    F_filt_blocked = [lowpass_filter(item, fs, cutoff=0.5) for item in F_blocked]

    # suite2p method
    dfof_method = 'suite2p'
    filter_first = False

    if dfof_method == 'suite2p':
        if filter_first:
            F_filt_blocked = [lowpass_filter(item, fs, cutoff=0.2) for item in F_blocked]

            dfof_s2p_ops = dict(baseline='constant_prctile', win_baseline=30, sig_baseline=5, prctile_baseline=20)
            blocked_fc = [dcnv.preprocess(F=item, fs=fs, **dfof_s2p_ops) for item in F_filt_blocked]

            blocked_zdfof = [zscore(item, axis=1) for item in blocked_fc]
            zdfof = np.hstack(tuple(blocked_zdfof)).astype('float32')

        else:
            dfof_s2p_ops = dict(baseline='maximin', win_baseline=45, sig_baseline=5, prctile_baseline=10)

            blocked_Fc = [dcnv.preprocess(F=item, fs=fs, **dfof_s2p_ops) for item in F_blocked]
            blocked_zdfof = [zscore(item, axis=1) for item in blocked_Fc]
            blocked_zdfof_filt = [lowpass_filter(item, fs, cutoff=0.5) for item in blocked_zdfof]

            F_corrected = np.hstack(tuple(blocked_Fc)).astype('float32')
            zdfof = np.hstack(tuple(blocked_zdfof)).astype('float32')
            zdfof_filt = np.hstack(tuple(blocked_zdfof_filt)).astype('float32')

#%%
    save_traces = True
    if save_traces:
        np.savez(f.with_name('RASTERPLOT_zdfof.npz'), F_corrected=F_corrected,
                 zdfof=zdfof, zdfof_filt=zdfof_filt
                 )


#%%

    # percentile method
    # elif dfof_method == 'prctile':
    #     dfof_prctile_ops = dict(win=90, percentile=10)
    #     blocked_dfof, blocked_df, blocked_baseline = compute_blocked_prctile_dfof(Ffilt_blocked, fs, **dfof_prctile_ops)
    #     blocked_zdfof = [zscore(item, axis=1) for item in blocked_dfof]
    #
    #     dfof = np.hstack(tuple(blocked_dfof)).astype('float32')
    #     zdfof = np.hstack(tuple(blocked_zdfof)).astype('float32')
    #
    # elif dfof_method == 'both':
    #     F_filt_blocked = [lowpass_filter(item, fs, cutoff=0.5) for item in F_blocked]
    #
    #     dfof_prctile_ops = dict(win=90, percentile=10)
    #     blocked_dfof, blocked_df, blocked_baseline = compute_blocked_prctile_dfof(F_filt_blocked, fs,
    #                                                                               **dfof_prctile_ops)
    #
    #     dfof_s2p_ops = dict(baseline='maximin', win_baseline=30, sig_baseline=5, prctile_baseline=20)
    #     blocked_fc = [dcnv.preprocess(F=item, fs=fs, **dfof_s2p_ops) for item in blocked_dfof]
    #
    #     blocked_zdfof = [zscore(item, axis=1) for item in blocked_fc]
    #     zdfof = np.hstack(tuple(blocked_zdfof)).astype('float32')

    # %%
    #
    # noise_levels = calculate_noise_levels(zdfof, fs)
    #
    # fig_noise, axarr = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    # sns.distplot(noise_levels, kde=True, ax=axarr[0])
    # axarr[0].set_title('noise level -  kde distplot')
    # sns.ecdfplot(noise_levels, ax=axarr[1])
    # axarr[1].set_title('noise level ecdf')
    # fig.suptitle(str(ainfo)
    #              + f"\nconc={blk_conc}\nreared={rear_type}"
    #              + f"\nnoise_level(dzfof)")
    # plt.xlim(0, 100)
    # plt.show()

    # # save processed traces, and noise figure
    # np.savez(f.with_name('zdfof.npy'), zdfof=zdfof, dfof_ops=dfof_ops, noise_levels=noise_levels)
    # fig_noise.savefig(f.with_name('fig_noise_zdfof.png'))
    # fig_noise.savefig(f.with_name('fig_noise_zdfof.pdf'))

    # %% align trials into 3d tensor (cells x trial x time
    zdfof_tensor, trial_ts = represp.align_signals(zdfof_filt, timestamps, df_stimulus.stim_ict)
    ustim = df_stimulus.stimname.unique().categories

    # calculate means across odors
    mean_odors = []
    for iodor, odor in enumerate(ustim):
        mean_odors.append(np.mean(zdfof_tensor[:, df_stimulus.stimname == odor, :], axis=1))

    mean_odors_zdfof = np.stack(mean_odors)
    mean_zdfof_hcat = np.hstack(mean_odors)

    # save trial tensors
    np.savez(f.with_name('RASTERPLOT_zdfof_tensor.npz'),
             zdfof_tensor=zdfof_tensor, trial_ts=trial_ts,
             mean_zdfof_tensor=mean_odors_zdfof,
             ustim=ustim.to_list())


    # %% SELECT ONLY GOOD CELLS TO RUN RASTERMAP
    mean_zdfof_hcat_rastermap = mean_zdfof_hcat[iscell, :]
    mean_odors_zdfof_raster = mean_odors_zdfof[:, iscell, :]
    # --------------------------------------------------------
    # run rastermap
    mean_model = Rastermap(**rasterplot_ops)
    mean_odor_embedding = mean_model.fit_transform(mean_zdfof_hcat_rastermap)

    # save mean odor embedding
    # --------------------------------------------------------
    np.save(f.with_name('embedding_mean_odor.npy'), dict(mean_zdfof_hcat=mean_zdfof_hcat_rastermap,
            mean_model=mean_model, mean_odor_embedding=mean_odor_embedding))

    isort1 = mean_model.isort
    stim_loc = np.argwhere(trial_ts == 0)[0, 0]

    fig_mean_raster, axarr = plt.subplots(1, len(ustim), figsize=(12, 12))
    for iodor, odor in enumerate(ustim):
        img = mean_odors_zdfof_raster[iodor, isort1, :]
        img_ax = axarr[iodor].imshow(img, cmap='gray_r', vmin=-0.5, vmax=1.5, aspect='auto')
        axarr[iodor].set_title(odor)
        axarr[iodor].get_xaxis().set_visible(False)
        axarr[iodor].get_yaxis().set_visible(False)
        axarr[iodor].axvline(stim_loc, linestyle='--', linewidth=1, color='r')
    axarr[0].get_yaxis().set_visible(True)
    fig_mean_raster.suptitle(
        str(ainfo) + f"\nconc={blk_conc}\nreared={rear_type}" + "\n" + 'averaged across odors, rastermap sorted')
    fig_mean_raster.colorbar(img_ax)
    plt.show()

    # save mean raster plots
    # --------------------------------------------------------
    fig_mean_raster.savefig(f.with_name('zdfof_odor_mean_rasterplot.pdf'))
    fig_mean_raster.savefig(f.with_name('zdfof_odor_mean_rasterplot.png'))



    # %%
    # cid0 = 20
    #
    # fig, axarr = vis.blocked_stackplot(zdfof_filt, timestamps, cids=np.arange(20) + cid0,
    #                                    df_stimulus=df_stimulus,
    #                                    df_frametimeinfo=df_frametimeinfo,
    #                                    peak_start=1, peak_len=5)
    #
    # if dfof_method == 'prctile':
    #     fig.suptitle(json.dumps(dfof_prctile_ops), fontsize=9)
    # elif dfof_method == 'suite2p':
    #     fig.suptitle(json.dumps(dfof_s2p_ops), fontsize=10)
    # elif dfof_method == 'both':
    #     fig.suptitle('both')
    # plt.show()

    # %% Rasterplot - full
    # -----------------------------------------------
    zdfof_filt_raster = zdfof_filt[iscell, :]
    model = Rastermap(**rasterplot_ops)
    embedding = model.fit_transform(zdfof_filt_raster)
    isort1 = np.argsort(embedding[:, 0])

    img = rasterplot.prep_to_plot(zdfof_filt_raster[isort1, :])
    vmax = np.max(img) * .8

    fig1, axarr1 = rasterplot.suite2p_display(img, isort=None, nbin=2, figure_kwargs={}, imshow_kwargs={})
    plt.show()

    fig2, axarr2 = rasterplot.suite2p_reordered_trial_display(img, df_stimulus, timestamps, nbin=2,
                                                              imshow_kwargs=dict(vmin=0, vmax=3))

    fig2.suptitle(title_str)
    plt.show()
    fig2.savefig(f.with_name('Rasterplot_reordered.pdf'))
    fig2.savefig(f.with_name('Rasterplot_reordered.png'))

    np.save(f.with_name('embedding_iscell.npy'), dict(zdfof_fill_iscell=zdfof_filt_raster,
            mean_model=model, embedding=embedding))
    # %%
    xid = model.xid
    uxid = np.unique(model.xid)

    zdfof_filt_raster

    clustered_responses = np.empty((uxid.size, zdfof_filt_raster.shape[1]))
    for i, clust in enumerate(np.unique(model.xid)):
        clustered_responses[i, :] = np.mean(zdfof_filt_raster[xid == clust, :], axis=0)

    img = rasterplot.prep_to_plot(clustered_responses)
    img = clustered_responses
    vmax = np.max(img) * .7
    fig, axarr = rasterplot.suite2p_reordered_trial_display(clustered_responses, df_stimulus, timestamps, nbin=1,
                                                            imshow_kwargs=dict(vmin=-vmax, vmax=vmax,
                                                                               interpolation=None,
                                                                               cmap='RdBu_r'))
    fig.suptitle(
        str(ainfo) + f"\nconc={blk_conc}\nreared={rear_type}" + "\n" + 'clustered responses from rastermap.')
    axarr[0].get_yaxis().set_visible(False)
    fig.savefig(f.with_name('rastermap_cluster_means.png'))
    fig.savefig(f.with_name('rastermap_cluster_means.pdf'))
    plt.show()

# %%
