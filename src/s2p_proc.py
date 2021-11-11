import itertools
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats, signal
from scipy.spatial.distance import pdist, squareform
from sklearn import preprocessing
from sklearn import svm
from dataload import SparsityAnalysis, load_all_flydata, make_fly_list, adjust_for_bad_frames
from rearing import hallem, io, corrmat, fti, tidy
from rearing.ca import respvec, trace, vis
from pandas.api.types import CategoricalDtype

import pyarrow.feather as feather
from suite2p.extraction import dcnv

from scipy.signal import convolve

Hallem = hallem.HallemSet()

all_stim_ord = ['pfo', '1-5ol', '1-6ol', 'EP', '1-6ol+EP', 'VA']

db = io.load_db()
db = db.loc[db['run']]
fly_list = make_fly_list(db)


def calc_movie_conc(df_stimulus):
    mconc = df_stimulus['conc'].to_numpy()
    mconc = np.unique(mconc[mconc < 0])[0].tolist()
    return mconc


def s2p_dfof(Fc):
    # baseline operation
    # get spikes

    df = dcnv.preprocess(
        F=Fc,
        baseline=ops['baseline'],
        win_baseline=60,
        sig_baseline=10,
        fs=ops['fs'],
        prctile_baseline=20
    )

    baseline = Fc - df
    dfof = df / baseline
    return dfof
    #
    # spks = dcnv.oasis(F=Fc, batch_size=ops['batch_size'], tau=ops['tau'], fs=ops['fs'])
    #
    # efilt = np.exp(- np.linspace(0, 50, 200) / (ops['tau'] * ops['fs']))
    # efilt /= efilt.sum()
    # sout = convolve(spks[0, :], efilt)
    # sout = sout[:spks.shape[1]]
    #
    # fig = px.line(dfof_s2p[:10, :].T, facet_row='variable')
    # fig.show(renderer='browser')


def process_experiment(mov, frametime_ops=None, filt_ops=None, dfof_ops=None, respvec_ops=None):
    ti = mov.load_ti()
    df_stimulus0 = mov.load_df_stimulus()
    df_frametimeinfo0 = mov.load_frametimeinfo()

    stim_list = df_frametimeinfo0.trialtrace_idx[df_frametimeinfo0.trialtrace_idx > 0].unique()

    s2p_dat0 = mov.load_suite2p_data()
    df_frametimeinfo, df_stimulus, s2p_dat = adjust_for_bad_frames(df_frametimeinfo0, df_stimulus0, s2p_dat0)

    if frametime_ops is not None:
        frame_times = df_frametimeinfo['frame_times'].to_numpy()
        df_frametimeinfo['peakresp_idx'] = fti.compute_peakresp_idx(frame_times,
                                                                    ti['stim_ict'],
                                                                    frametime_ops['peakresp_win'])

    # --------------------------------------------------------------
    # DF_STIMULUS - make 'stimname' a categorical variable
    # --------------------------------------------------------------
    stim_ord = [item for item in all_stim_ord if item in df_stimulus['stimname'].unique()]  # stimname ordering
    cat_type = CategoricalDtype(categories=stim_ord, ordered=True)

    df_stimulus['stimname'] = pd.Categorical(df_stimulus['stimname'], ordered=True, categories=stim_ord)

    # --------------------------------------------------------------
    # DF_STIMULUS - add scope_pulse_idx, stim_ict, stim_fct
    # --------------------------------------------------------------
    df_stimulus['block_idx'] = df_frametimeinfo.loc[df_frametimeinfo['olf_pulse_idx'] > 0, :].groupby('olf_pulse_idx') \
        .first()['scope_pulse_idx'].to_list()

    df_stimulus['stim_ict'] = ti['stim_ict']
    df_stimulus['stim_fct'] = ti['stim_fct']
    # --------------------------------------------------------------
    # compute movie concentration
    # --------------------------------------------------------------
    movie_conc = calc_movie_conc(df_stimulus)

    # --------------------------------------------------------------
    # compute frame rate
    # --------------------------------------------------------------
    fs = trace.compute_fs(df_frametimeinfo['frame_times'].to_numpy())

    # --------------------------------------------------------------
    # default parameters
    # --------------------------------------------------------------
    if filt_ops is None:
        filt_ops = dict(fs=fs, cutoff=0.5, trans_width=0.05, numtaps=230)
    if dfof_ops is None:
        dfof_ops = dict(win=45, percentile=20)

    params = dict(frametime_ops=frametime_ops,
                  filt_ops=filt_ops,
                  dfof_ops=dfof_ops,
                  stim_ord=stim_ord,
                  movie_conc=movie_conc)

    trials_and_timing = dict(df_frametimeinfo=df_frametimeinfo,
                             df_stimulus=df_stimulus,
                             ti=ti)

    # --------------------------------------------------------------
    # low pass filter traces
    # --------------------------------------------------------------
    blocked_traces = trace.split_traces_to_blocks(s2p_dat['Fc'], df_frametimeinfo.scope_pulse_idx)
    b = signal.remez(filt_ops['numtaps'], [0, filt_ops['cutoff'], filt_ops['cutoff'] + filt_ops['trans_width'],
                                           0.5 * filt_ops['fs']], [1, 0], Hz=filt_ops['fs'])
    filt_blocked_traces = [signal.filtfilt(b, 1, y) for y in blocked_traces]

    # --------------------------------------------------------------
    # extract dF/F
    # --------------------------------------------------------------
    dfof_method = "suite2p"
    if dfof_method == "remy":
        blocked_dfof, blocked_df, blocked_baseline = zip(*[trace.extract_dfof(item,
                                                                              fs=fs,
                                                                              win=dfof_ops['win'],
                                                                              percentile=dfof_ops['percentile'])
                                                           for item in filt_blocked_traces])
    elif dfof_method == "suite2p":
        blocked_dfof, blocked_df, blocked_baseline = zip(*[trace.suite2p_extract_dfof(
            F=item,
            baseline=s2p_dat['ops']['baseline'],
            win_baseline=60,
            sig_baseline=10,
            fs=s2p_dat['ops']['fs'],
            prctile_baseline=10) for item in blocked_traces])


    dfof = np.hstack(blocked_dfof).astype('float32')
    df = np.hstack(blocked_df).astype('float32')
    baseline = np.hstack(blocked_baseline).astype('float32')
    traces = dict(dfof=dfof, df=df, baseline=baseline)
    # --------------------------------------------------------------
    # compute response vectors
    # --------------------------------------------------------------
    # baseline avg & mean, peakresp avg, peakresp amplitude
    peakresp_nmean, baseline_avg, baseline_std = \
        respvec.compute_response_vectors(dfof,
                                         stim_list=stim_list,
                                         baseline_idx=df_frametimeinfo.baseline_idx.to_numpy(),
                                         peakresp_idx=df_frametimeinfo.peakresp_idx.to_numpy(),
                                         maxk=respvec_ops['maxk'],
                                         stack_numpy=True)
    peakresp_amp = peakresp_nmean - baseline_avg
    dist, df_peakresp_amp = respvec.compute_corrmats(peakresp_amp, df_stimulus=df_stimulus)

    # calculate peak locations
    peakresp_max, peakmax_idx = respvec.find_peaks_and_idx(dfof,
                                                           df_frametimeinfo.peakresp_idx.to_numpy(),
                                                           stack_numpy=True)
    peakmax_times = np.vstack([frame_times[row] for row in peakmax_idx])

    respvecs = dict(peakresp_amp=peakresp_amp,
                    peakresp_nmean=peakresp_nmean,
                    baseline_avg=baseline_avg,
                    baseline_std=baseline_std,
                    peakresp_max=peakresp_max,
                    peakmax_idx=peakmax_idx,
                    peakmax_times=peakmax_times,
                    df_peakresp_amp=df_peakresp_amp,
                    dist=dist
                    )

    if respvec_ops is not None:
        bin_response = respvec.get_binresponse(peakresp_amp, baseline_std,
                                               n_std=respvec_ops['n_std'], thr=respvec_ops['thr'])
        df_tuning, df_binresp = respvec.get_tuning(bin_response, df_stimulus)
        df_tuning_logical = (df_tuning > 0) * 1
        binresp = dict(bin_response=bin_response,
                       df_binresp=df_binresp,
                       df_tuning=df_tuning,
                       df_tuning_logical=df_tuning_logical,
                       ops=respvec_ops)
        params['respvec_ops'] = respvec_ops
    else:
        binresp = dict()

    return params, trials_and_timing, traces, respvecs, binresp


def param_to_subdir(params):
    subdirnam = f"__win{params['dfof_ops']['win']}__prctile{params['dfof_ops']['percentile']}" \
                f"__peakwin{params['frametime_ops']['peakresp_win']}"
    return subdirnam

    # def plot_dfof():

    """"""

    # cids = list(range(n_cells))
    #
    # tidy_resppeaks = tidy.tidy_response_peaks(respvecs['peakresp_max'],
    #                                           respvecs['peakmax_times'],
    #                                           bin_response=,
    #                                           df_frametimeinfo=df_frametimeinfo)
    # #
    # with sns.axes_style('white', {"xtick.major.size": 5, "ytick.major.size": 5, "axes.linewidth": 0.5,
    #                               "axes.edgecolor": '0.8'}):
    #     cid0 = 0
    #     plot_cids = list(range(cid0, cid0 + 20))
    #
    #     g = sns.FacetGrid(data=tidy_dfof[tidy_dfof['cid'].isin(plot_cids)], row='cid', col='scope_pulse_idx',
    #                       ylim=[0, 2],
    #                       margin_titles=True,
    #                       sharex='col', sharey=True, aspect=9, height=0.5)
    #     g.map_dataframe(sns.lineplot, x='frame_times', y='dfof', clip_on=False, zorder=2.5)
    #     g.set_titles(col_template="block {col_name}", row_template="{row_name}")
    #     g.fig.subplots_adjust(top=0.9, right=0.95)
    #     g.fig.suptitle('df/f traces')
    #
    #     # mark responders
    #     for (row_val, col_val), ax in g.axes_dict.items():
    #         mask = (tidy_resppeaks['cid'] == row_val) & (tidy_resppeaks['scope_pulse_idx'] == col_val)
    #         df = tidy_resppeaks.loc[mask, :]
    #         ax.plot(df.peakmax_time, df.peakresp_max + 0.5, 'rv',
    #                 markersize=3, clip_on=False)
    #         ax.set_axisbelow(True)
    #
    #     # shade stimulus intervals
    #     for i, ax in enumerate(g.axes[0, :]):
    #         df_shade = df_stimulus.loc[df_stimulus['block_idx'] == i + 1, :]
    #         vis.draw_stimuli(ax, df_shade['stim_ict'], df_shade['stim_fct'], df_shade['stimname'])
    #
    #     sns.despine(offset=5, trim=False)
    #     trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    #     plt.show()
    return None


# def main(fly_num, expt_num):
#     fly = fly_list[fly_num]
#     mov = fly.expt_list[expt_num]
#
#     return fly, mov, params, trials_and_timing, traces, respvecs


def multiconcentration_corrmats(movies, frametime_ops, respvec_ops, dfof_ops):
    """

    Parameters
    ----------
    movies : list of experiments to include

    frametime_ops : dict of parameters for dfof extraction and peak amplitude quantification
                    example - {'peakresp_win': 15, 'baseline_win': -15, 'trial_win': 20}

    respvec_ops : dict of parameters for responder criteria
                    example - {'n_std': 5, 'thr': 0.5, 'maxk': 20}

    Returns
    -------
        fig - plotly chart
fig = multiconcentration_corrmats(mov, frametime_ops, respvec_ops)
    """
    correlogram_list = []
    for mov in movies:
        params, trials_and_timing, traces, respvecs, binresp = process_experiment(mov,
                                                                                  frametime_ops=frametime_ops,
                                                                                  respvec_ops=respvec_ops,
                                                                                  dfof_ops=dfof_ops)

        tstr = f"date: {mov.fly_parent.date}, fly: {mov.fly_parent.fly_num}, movies: {mov.movie} [{params['movie_conc']}] " \
               f"\nrearing: {mov.fly_parent.metadata['Fly info'][mov.fly_parent.fly_num]['reared']}"

        fig = vis.correlograms(respvecs['dist'], title=tstr)

        correlogram_list.append(fig)
    return correlogram_list




def corrmat_violins(movies, frametime_ops, respvec_ops):
    tidy_tril_list = []
    for mov in movies:
        params, trials_and_timing, traces, respvecs, binresp = process_experiment(mov,
                                                                                  frametime_ops=frametime_ops,
                                                                                  respvec_ops=respvec_ops)
        tidy_tril = tidy.tidy_tril_df_corrmat(respvecs['dist']['pearson'])
        tidy_tril['movie_conc'] = params['movie_conc']
        tidy_tril_list.append(tidy_tril)
    return pd.concat(tidy_tril_list)


if __name__ == '__main__':

    figures_to_make = ['corrmats']

    for fly in fly_list:
        for mov in fly.expt_list:
            # options
            frametime_ops = dict(peakresp_win=15, baseline_win=-15, trial_win=20)
            respvec_ops = dict(n_std=5, thr=0.5, maxk=10)
            dfof_ops = dict(win=60, percentile=20)

            params, trials_and_timing, traces, respvecs, binresp = process_experiment(mov,
                                                                                      frametime_ops=frametime_ops,
                                                                                      respvec_ops=respvec_ops,
                                                                                      dfof_ops=dfof_ops)

            vis.corrmats(respvecs['dist'])
            plt.show()

            tidy_tril = tidy.tidy_tril_df_corrmat(respvecs['dist']['pearson'])
            params['respvec_ops'] = respvec_ops

            dfof = traces['dfof']
            df_frametimeinfo = trials_and_timing['df_frametimeinfo']
            df_stimulus = trials_and_timing['df_stimulus']
            ti = trials_and_timing['ti']

            n_cells, n_ts = dfof.shape
            n_trials = df_stimulus.shape[0]
            params['n_cells'] = n_cells

            tidy_dfof = tidy.tidy_traces(dfof, df_frametimeinfo)

            # ---------------------------------------------------------------------------
            # save analysis outputs
            # ---------------------------------------------------------------------------
            savesubdir = mov.PROC_DIR().joinpath(param_to_subdir(params))
            if not savesubdir.exists():
                savesubdir.mkdir()

            # save params to yaml file
            with open(savesubdir.joinpath('params.yaml'), 'w') as outfile:
                yaml.dump(params, outfile, default_flow_style=False, sort_keys=False)

            # save df_stimulus
            df_stimulus.to_csv(savesubdir.joinpath('df_stimulus.csv'), index=False)

            # save .npy files
            np.save(savesubdir.joinpath('params.npy'), params)
            np.save(savesubdir.joinpath('trials_and_timing.npy'), trials_and_timing)
            np.save(savesubdir.joinpath('traces.npy'), traces)
            np.save(savesubdir.joinpath('respvecs.npy'), respvecs)
            if len(binresp) > 0:
                np.save(savesubdir.joinpath('binresp.npy'), binresp)

            feather.write_feather(tidy_dfof, savesubdir.joinpath('tidy_dfof.feather'), compression='uncompressed')
            print('\n==============================')
            print(mov)
            print('.........saved files.')
            # ---------------------------------------------------------------------------
            # corrmats
            # ---------------------------------------------------------------------------
            tstr = f"date: {fly.date}, fly: {fly.fly_num}, movies: {mov.movie} [{params['movie_conc']}] " \
                   f"\nrearing: {fly.metadata['Fly info'][fly.fly_num]['reared']}"

            if figures_to_make is not None:
                if 'corrmats' in figures_to_make:
                    fig_corrmat, axarr_corrmat, haxarr_corrmat = vis.corrmats(respvecs['dist'])
                    fig_corrmat.suptitle(tstr)
                    fig_corrmat.savefig(savesubdir.joinpath('corrmats.pdf'))
                    fig_corrmat.savefig(savesubdir.joinpath('corrmats.png'))

        # ------------------------------------------------------------------
        # plot traces, with responders and stimuli showngi
        # ------------------------------------------------------------------

        print('done')


