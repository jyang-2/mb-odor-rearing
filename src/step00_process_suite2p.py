"""
Run after everything in preprocess_suite2p.py is complete.

Saves C_dec to the same folder as the suite2p statfile.

"""
import copy
import itertools
import json
# from run_correlation_analysis import mean_corrmats
import pprint as pp
from collections import namedtuple

import caiman
import caiman.source_extraction.cnmf as cm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn import preprocessing
from suite2p.extraction import dcnv

from rearing import tidy, hallem
from rearing.ca import respvec, trace, represp
from rearing_datasets import get_pfo_statfiles, get_odor_statfiles, get_pfo_datasets, get_odor_datasets
from ryim import dirutils, analysis

mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 10, "savefig.format": 'png', 'text.usetex': False})



#  Set ordering of odors
#  Use to create a categorical dtype (ordered) in pandas. This can be used to set the dtype in df_stimulus['stimname'],
#  for easy ordering.
#  To order by stimname, do df_stimulus.sort_values(by='stimname')
odor_cattype = pd.CategoricalDtype(categories=['pfo', '1-5ol', 'VA', '1-6ol', 'EP', '1-6ol+EP'], ordered=True)

# namedtuple, structure for holding metadata
MetaData = namedtuple("MetaData",
                      ["df_stimulus",
                       "df_frametimeinfo",
                       "timestamps",
                       "session_info",
                       "xml_meta",
                       "sync_meta",
                       "meta_yaml",
                       "fs", "blk_conc", "rear_type"
                       ])

Hallem = hallem.HallemSet()

readme_txt = ("\n"
              "This folder (step00_process_suite2p) contains the outputs of running step00_process_suite2p.py on the initial suite2p outputs.\n"
              "\n"
              "It should contain the following files:\n"
              "	- params.json: contains the parameters used for extracting C_dec, the cells used, and the cost\n"
              "			- dcnv_params\n"
              "			- oasis_params\n"
              "			- respvec_ops \n"
              "			- good_cells: cells used to compute cost function, specified by user\n"
              "	\n"
              "	- model_selection.json: contains a list of \"models\" (or parameter combos) to try, and their cost. \n"
              "							The selected parameters in params.json come from the \"model\" with the lowest cost.\n"
              "	\n"
              "	- Fcc: outputs from baseline correction w/ suite2p.dcnv, run block by block (using df_frametimeinfo) before np.hstack\n"
              " - Fq: quantile-normed Fcc (preprocess.RobustScaler)\n"
              " - C_dec: deconvolved traces\n"
              " - Sp : Spike array"
              "	- oasis_results: final output of deconvolution\n"
              "	-    \n")


# %%
def correct_and_deconvolve(Fc, df_frametimeinfo, df_stimulus, dcnv_params, oasis_params):
    """
        Takes parameters for baseline correction (suite2p.dcnv.preprocess) and deconvolution (oasis AR2, constrained
        foopsi) and returns block-corrected traces.

        :param Fc: neuropil-corrected fluorescence from suite2p ( = F - 0.7 * Fneu)
        :type Fc: np.ndarray
    """
    n_cells, T = Fc.shape
    Fc = trace.split_traces_to_blocks(Fc, df_frametimeinfo.scope_pulse_idx)

    Fcc = [dcnv.preprocess(F=item, **dcnv_params) for item in Fc]
    Fcc = np.hstack(tuple(Fcc))

    # quantile matching & transform, quantile_range=(0.25, 0.75)
    d = preprocessing.RobustScaler().fit_transform(Fcc)

    # oasis deconvolution
    C_dec = np.zeros(d.shape)
    Bl = np.zeros(n_cells)
    C1 = np.zeros(d.shape)
    G = np.zeros((n_cells, 2))
    Sn = np.zeros(n_cells)
    Sp = np.zeros((n_cells, T))
    Lam = np.zeros(n_cells)

    for cid in range(C_dec.shape[0]):
        y = d[cid, :]

        g = caiman.source_extraction.cnmf.deconvolution.estimate_time_constant(y, p=2, lags=5, fudge_factor=1.)
        c, bl, c1, g, sn, sp, lam = caiman.source_extraction.cnmf.deconvolution.constrained_foopsi(y, **oasis_params)
        C_dec[cid, :] = c
        Bl[cid] = bl
        C1[cid] = c1
        G[cid, :] = g
        Sp[cid, :] = sp
        Sn[cid] = cm.deconvolution.GetSn(y)
        Lam[cid] = lam

    # save oasis results
    oasis_results = dict(C_dec=C_dec, bl=Bl, c1=C1, g=G, sp=Sp, lam=Lam, oasis_params=oasis_params)

    return Fcc, d, oasis_results


def cost_fcn(Fc_, timestamps_, df_frametimeinfo_, df_stimulus_, dcnv_params_, oasis_params_, respvec_ops_):
    """
    Function used to pick the parameters for processing fluorescence traces.

    :param Fc_: 2d np.array, [cells x time], F - 0.7*Fneu
    :param timestamps_: 1d np.array
    :param df_frametimeinfo_:
    :param df_stimulus_:
    :param dcnv_params_: parameters for suite2p.dcnv.preprocess
    :param oasis_params_: parameters for caiman.source_extraction.cnmf.deconvolution.constrained_foopsi
    :param respvec_ops_: parameters for calculating response vectors
    :return:
        - cost: geometric mean of standard deviation of intra-odor response correlations
    """
    Fcc_, d_, oasis_results_ = correct_and_deconvolve(Fc_, df_frametimeinfo_, df_stimulus_, dcnv_params_, oasis_params_)

    C_dec_ = oasis_results_['C_dec']

    # compute response magnitudes
    peak_mean_, baseline_mean_, baseline_std_ = respvec.compute_interval_response_vectors(C_dec_,
                                                                                          timestamps_,
                                                                                          df_stimulus_['stim_ict'],
                                                                                          **respvec_ops_)

    peak_amp_ = peak_mean_ - baseline_mean_

    dist_all_, df_peakresp_ = respvec.compute_corrmats(peak_amp_, df_stimulus_)
    # dist_all_mean_ = mean_corrmats(dist_all_)

    tidy_corr_ = tidy.tidy_tril_df_corrmat(dist_all_['pearson']).drop(['stim_idx_col', 'stim_idx_row'], axis=1)
    grouped_ = tidy_corr_.groupby(by=['stimname_row', 'stimname_col'])
    cost_ = stats.gmean(grouped_.std().dropna().to_numpy())[0]
    return cost_


def get_title_str(ainfo_, nt_meta_data_):
    """
    Constructs title for plotting/analysis

    :param ainfo_: info returned by dirutils.parse_dataset_info('.../suite2p/**/stat.npy')
    :type ainfo_: namedtuple AnalysisSet in ryim.dirutils
    :param nt_meta_data_: holds timing, stimulus, recording info
    :type nt_meta_data_: namedtuple MetaData
    :return: f-string, with relevant movie/block information
    :rtype: str

    """

    title_str = f"""{ainfo_.flydate} fly {ainfo_.flynum}, {ainfo_.mov}\n"""

    if ainfo_.is_block:
        title_str = f"""{title_str} blk={ainfo_.blk}\n"""

    title_str = title_str + f"""reared: {nt_meta_data_.rear_type}    |   conc={nt_meta_data_.blk_conc}"""
    return title_str


def postprocess_metadata(_nt_meta_data: namedtuple, ocat=odor_cattype):
    """
    Fix some of the fields in _nt_meta_data.

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
    _df_stimulus['stimname'] = _df_stimulus['stimname'].astype(odor_cattype)

    _df_stimulus['reps'] = ''
    for ustim in _df_stimulus['stimname'].unique():
        mask = _df_stimulus['stimname'] == ustim
        _df_stimulus.loc[mask, 'reps'] = np.arange(np.count_nonzero(mask * 1))

    _session_info = {key: value for key, value in _nt_meta_data.session_info.to_dict().items() if key.isidentifier()}

    _fs = trace.compute_fs(_nt_meta_data.timestamps)
    _blk_conc = _df_stimulus.blk_conc.unique()[0]
    _rear_type = _session_info['rearing']

    return MetaData(_df_stimulus,
                    _df_frametimeinfo,
                    _nt_meta_data.timestamps,
                    _session_info,
                    _nt_meta_data.xml_meta,
                    _nt_meta_data.sync_meta,
                    _nt_meta_data.meta_yaml,
                    _fs,
                    _blk_conc,
                    _rear_type
                    )


def make_label_for_svm(ts, stim_ict, stim_idx, peak_start, peak_len):
    """
    Fill a vector (like timestamps) w/ stim_idx at specified intervals w/ after stim_ict.

    :param ts: timestamps
    :type ts: np.ndarray, 1D
    :param stim_ict: time that stimulus turns on; df_stimulus.stim_ict
    :type stim_ict: np.ndarray
    :param stim_idx: 1-indexed stimulus # index (ex: df_stimulus.stim_idx)
    :type stim_idx:  np.ndarray
    :param peak_start: time (in seconds) after stim_ict the peak window starts
    :type peak_start: float
    :param peak_len: time (in seconds) after stim_ict + peak_start, how long the peak window lasts for
    :type peak_len: float
    :return: labeltimestamp
    :rtype:  np.array

    """
    label_ts = np.zeros_like(ts)

    for s_ict, s_idx in zip(stim_ict, stim_idx):
        t0 = s_ict + peak_start
        t1 = s_ict + peak_start + peak_len
        idx = (ts >= t0) & (ts < t1)
        label_ts[idx] = s_idx
    return label_ts


def draw_oasis_results(y_, c_, sp_, x_=None, sn_=None, ax_=None):
    """
    Plots fluorescence activity and results of calcium deconvolution.

    :param y_: fluorescence trace, raw
    :param c_: deconvolved trace
    :param sp_: spikes
    :param x_: timestamps for y_, c_, sp_
    :param sn_: SNR, scalar
    :param ax_: axes
    :return: axes w/ plots
    """
    lw = 0.5

    if x_ is None:
        x_ = np.arange(y_.size)

    if ax_ is None:
        ax_ = sns.lineplot(x=x_, y=y_, color='k', alpha=0.3, lw=0.5)
    else:
        sns.lineplot(x=x_, y=y_, color='k', alpha=0.3, ax=ax_, lw=0.5)

    sns.lineplot(x=x_, y=c_, color='k', alpha=1, ax=ax_, lw=0.5)
    sns.rugplot(x=x_[sp_ > 0], height=0.1, color='r', ax=ax_)

    if sn_ is not None:
        ax_.axhline(sn_, linestyle='--', alpha=0.5)

    return ax_


def load_metadata(f, method='statfile'):
    """Helper func"""
    if method == 'statfile':
        ainfo = dirutils.parse_dataset_info(f)

    elif method == 'AnalysisSet':
        ainfo = f

    print(ainfo)
    # loading metadata
    if ainfo.is_block:
        meta_data_loader = analysis.BlockMetaDataLoader.init_from_statfile(f)
    else:
        meta_data_loader = analysis.MovieMetaDataLoader.init_from_statfile(f)
    return meta_data_loader.load()


def load_suite2p_data(f, method='statfile'):
    """
    Helper func for eaasy loading of data

    :param f: Path to stat.npy
    :param method: {'statfile' or 'AnalysisSet'}

    """
    if method == 'statfile':
        ainfo = dirutils.parse_dataset_info(f)
    elif method == 'AnalysisSet':
        ainfo = f

    # loading suite2p output data
    if f.parent.name == 'combined':
        suite2p_loader = analysis.Suite2pCombinedLoader(statfile=f)
    elif 'plane' in f.parent.name:
        suite2p_loader = analysis.Suite2pPlane0Loader(statfile=f)

    nt_s2p_data_ = suite2p_loader.load()
    return nt_s2p_data_


def load_all(f):
    """
    Load metadata and suite2p outputs from movie/movie block for stat.npy file.

    :param f: Path to stat.npy
    :return:
        nt_meta_data
        nt_s2p_data:
        title_str:
    """
    ainfo = dirutils.parse_dataset_info(f)
    nt_meta_data = postprocess_metadata(load_metadata(f))
    nt_s2p_data = load_suite2p_data(f)
    title_str = get_title_str(ainfo, nt_meta_data)
    return nt_meta_data, nt_s2p_data, title_str


def main(f):
    ainfo = dirutils.parse_dataset_info(f)

    # get data and info for statfile
    nt_meta_data, nt_s2p_data, title_str = load_all(f)

    # %% extract loaded files
    """
    STEP 02: process calcium responses, make them into population vectors
    """
    df_stimulus = copy.deepcopy(nt_meta_data.df_stimulus)
    df_frametimeinfo = copy.deepcopy(nt_meta_data.df_frametimeinfo)
    timestamps = copy.deepcopy(nt_meta_data.timestamps)
    ti = copy.deepcopy(nt_meta_data.sync_meta)
    iscell = nt_s2p_data.iscell[:, 0].astype(np.bool_)
    cellprob = nt_s2p_data.iscell[:, 1]

    if nt_meta_data.xml_meta['z_steps'] == 1:
        use_saved_cells = False
    else:
        use_saved_cells = False

    if use_saved_cells:
        good_cells_bool = iscell
    else:
        good_cells_bool = (cellprob > np.quantile(cellprob, .5)).squeeze()
    good_cells = np.argwhere(good_cells_bool).squeeze()

    # %% compute quantile-corrected cell traces, cell ranking by probability
    Fc = nt_s2p_data.F - 0.7 * nt_s2p_data.Fneu
    n_cells, T = Fc.shape
    cell_ranking = np.argsort(cellprob)[::-1]

    # parameters for trace processing
    models = []

    for i, (win, prc, sig) in enumerate(itertools.product([45, 60, 90], [10], [3, 5])):
        print(f"{i} ...win={win}, prc={prc}")
        print(f"-----------------------------------")

        # iterate over parameters to try
        dcnv_params = dict(baseline='maximin',
                           win_baseline=win,
                           sig_baseline=sig,
                           fs=nt_meta_data.fs,
                           prctile_baseline=10)

        oasis_params = dict(p=2, penalty=0, g=np.array([1.1, -0.15]), s_min=1.0, g_optimize=3)

        respvec_ops = dict(peak_start=0.5, peak_len=5, baseline_len=-10)

        params = dict(dcnv_params=dcnv_params, oasis_params=oasis_params, respvec_ops=respvec_ops)
        pp.pprint(params)

        cost = cost_fcn(Fc[good_cells, :], nt_meta_data.timestamps, nt_meta_data.df_frametimeinfo,
                        nt_meta_data.df_stimulus,
                        dcnv_params, oasis_params, respvec_ops)

        print(f'cost = {cost}')

        models.append(dict(params=params, cost=cost))

    # serialize models for json
    for item in models:
        item['params']['oasis_params']['g'] = item['params']['oasis_params']['g'].tolist()

    # choose best parameters
    for i, m in enumerate(models):
        print(f"model {i}: cost={m['cost']}")
        print("\n")

    cost = np.array([item['cost'] for item in models])
    model_idx = np.argmin(cost)

    best_params = copy.deepcopy(models[model_idx]['params'])
    best_params['good_cells'] = good_cells.tolist()  # cost function calculated only for good cells

    # create folder in the same parent folder that holds stat.npy
    outputdir = f.parent.joinpath("step00_process_suite2p")
    outputdir.mkdir(parents=True, exist_ok=True)

    # save readme, for future self
    with open('readme.txt', 'w') as outfile:
        outfile.write(readme_txt)

    # save best parameters to params.json
    with open(outputdir.joinpath('params.json'), 'w') as outfile:
        json.dump(best_params, outfile, indent=4)

    # save results of model/parameter selection to param_selection.json
    with open(outputdir.joinpath('param_selection.json'), 'w') as outfile:
        json.dump(models, outfile, indent=4)

    # Run trace correction again, but on all cells, in case good_cells changes
    Fcc, Fq, oasis_results = correct_and_deconvolve(Fc, nt_meta_data.df_frametimeinfo, df_stimulus,
                                                    best_params['dcnv_params'],
                                                    best_params['oasis_params'])
    C_dec = oasis_results['C_dec']
    Sp = oasis_results['sp']

    # Save outputs of process
    np.save(outputdir.joinpath('Fcc.npy'), Fcc)
    np.save(outputdir.joinpath('Fq.npy'), Fq)
    np.save(outputdir.joinpath('C_dec.npy'), C_dec)
    np.save(outputdir.joinpath('Sp.npy'), Sp)
    np.save(outputdir.joinpath('oasis_results.npy'), oasis_results, allow_pickle=True)
    return True


if __name__ == '__main__':
    plt.close('all')

    pfo_datasets = get_pfo_datasets()
    odor_datasets = get_odor_datasets()

    pfo_stat_files = get_pfo_statfiles()
    odor_stat_files = get_odor_statfiles()

    for f in pfo_stat_files:
        main(f)

    for f in odor_stat_files:
        main(f)

# #%% correlation analysis, for sanity check
#
# # compute response magnitudes
# respvec_ops = dict(peak_start=0.5, peak_len=5, baseline_len=-10)
# peak_mean, baseline_mean, baseline_std = respvec.compute_interval_response_vectors(C_dec,
#                                                                                    timestamps,
#                                                                                    df_stimulus['stim_ict'],
#                                                                                    **respvec_ops)
#
# peak_amp = peak_mean - baseline_mean
# # %%
# dist_all, df_peakresp = respvec.compute_corrmats(peak_amp[good_cells, :], df_stimulus)
# dist_all_mean = mean_corrmats(dist_all)
#
# tidy_corr = tidy.tidy_tril_df_corrmat(dist_all['pearson']).drop(['stim_idx_col', 'stim_idx_row'], axis=1)
# grouped = tidy_corr.groupby(by=['stimname_row', 'stimname_col'])
# tidy_cost = grouped.agg(['mean', 'std']).dropna()
# cost = stats.gmean(grouped.std().dropna().to_numpy())
#
# fig_corrmat_all, axarr, haxarr = vis.corrmats(dist_all)
# fig_corrmat_all.suptitle(title_str)
# plt.show(block=False)
#
# fig_corrmat_respmean, axarr, haxarr = vis.mean_corrmats(dist_all_mean, annot=True)
# fig_corrmat_respmean.suptitle(title_str)
# plt.show(block=False)
# # %%
# def plot_cell(cid_):
#     ax = draw_oasis_results(d[cid_, :], C_dec[cid_, :], Sp[cid_, :],
#                             sn_=Sn[cid_], x_=timestamps)
#     ax.set_title(f"cell # {cid_} (snr={round(Sn[cid_], 2)})", fontsize=8)
#     ax.spines[["left", "bottom"]].set_position(("data", 0))
#     plt.subplots_adjust(bottom=0.05, top=0.9, hspace=0.25)
#     plt.suptitle(title_str, fontsize=10)
#     plt.show(block=False)
#
#
# # %% plot oasis results
# save_oasis_plots = True
# if save_oasis_plots:
#     # of cells per figure
#     n_axes = 40
#     for cid0 in range(0, n_axes, 41):
#         fig, axarr = plt.subplots(10, int(np.ceil(n_axes / 10)), figsize=(11, 8),
#                                   sharey='all', sharex='all')
#
#         for cid, ax in zip(range(cid0, cid0 + n_axes), axarr.flatten()):
#             ax = draw_oasis_results(d[cid, :], C_dec[cid, :], Sp[cid, :],
#                                     sn_=Sn[cid], x_=timestamps, ax_=ax)
#
#             ax.set_title(f"cell # {cid} (snr={round(Sn[cid], 2)})", fontsize=8)
#             ax.spines[["left", "bottom"]].set_position(("data", 0))
#             # ax.set_ylim([-2, 10])
#
#     plt.subplots_adjust(bottom=0.05, top=0.9, hspace=0.25)
#     plt.suptitle(title_str, fontsize=10)
#     plt.show(block=False)
