import numpy as np
import heapq
import pandas as pd
from scipy.spatial.distance import pdist, squareform


def get_tuning(bin_response, df_stimulus):
    df_binresp = pd.DataFrame(bin_response)
    df_binresp.columns = pd.MultiIndex.from_frame(df_stimulus.loc[:, ['stim_idx', 'stimname']])
    df_tuning = df_binresp.drop(columns=[(1, 'pfo')]).groupby(by='stimname', axis=1).sum()
    return df_tuning, df_binresp


def get_binresponse(peakresp_amp, baseline_std, n_std=3, thr=0.5):
    binresp = (peakresp_amp >= baseline_std * n_std) & (peakresp_amp >= thr)
    return binresp


def compute_response_baselines(dfof, baseline_idx, stim_list=None, stack_numpy=True, verbose=False):
    """
    :param stim_list: 1D numpy array (or list) of stim. #s to use for indexing into baseline_idx
    :type stim_list: np.array
    """
    if stim_list is None:
        stim_list = np.unique(baseline_idx).nonzero()[0]

    n_trials = stim_list.size
    baseline_std = [None] * n_trials
    baseline_avg = [None] * n_trials
    for i, stim_idx in enumerate(stim_list):
        idx = (baseline_idx == stim_idx)
        # baseline mean
        bavg = np.mean(dfof[:, idx], axis=1)
        baseline_avg[i] = bavg
        # baseline stdev
        bstd = np.std(dfof[:, idx], axis=1)
        baseline_std[i] = bstd

    if stack_numpy:  # stack list of numpy arrays into 1 array
        baseline_std = np.stack(baseline_std, axis=1)
        baseline_avg = np.stack(baseline_avg, axis=1)

    return baseline_avg, baseline_std


def compute_response_peaks(dfof, peakresp_idx, maxk=3, stim_list=None, stack_numpy=True, verbose=False):
    if stim_list is None:
        stim_list = np.unique(peakresp_idx[peakresp_idx > 0])
        # stim_list = np.unique(peakresp_idx).nonzero()[0]

    n_trials = stim_list.size
    peakresp_nmean = [None] * n_trials

    for i, stim_idx in enumerate(stim_list):
        idx = np.nonzero(peakresp_idx == stim_idx)[0]
        peakresp = dfof[:, idx]
        npk = np.array([np.mean(heapq.nlargest(maxk, row)) for row in peakresp])
        peakresp_nmean[i] = npk
    if stack_numpy:
        peakresp_nmean = np.stack(peakresp_nmean, axis=1)
    return peakresp_nmean


def find_peaks_and_idx(dfof, peakresp_idx, stim_list=None, stack_numpy=True):
    if stim_list is None:
        stim_list = np.unique(peakresp_idx[peakresp_idx > 0])

    n_trials = stim_list.size
    peakresp_max = [None] * n_trials
    peakmax_idx = [None] * n_trials
    for i, stim_idx in enumerate(stim_list):
        idx = np.nonzero(peakresp_idx == stim_idx)[0]
        peakresp = dfof[:, idx]
        maxpeak = np.amax(peakresp, axis=1)
        maxidx = np.argmax(peakresp, axis=1)
        peakresp_max[i] = maxpeak
        peakmax_idx[i] = idx[maxidx]
    if stack_numpy:
        peakresp_max = np.stack(peakresp_max, axis=1)
        peakmax_idx = np.stack(peakmax_idx, axis=1)
    return peakresp_max, peakmax_idx


def compute_response_vectors(dfof, baseline_idx, peakresp_idx, stim_list=None, maxk=3, stack_numpy=True):
    if stim_list is None:
        stim_list = np.unique(baseline_idx).nonzero()[0]
    n_trials = stim_list.size

    baseline_avg, baseline_std = compute_response_baselines(dfof, baseline_idx=baseline_idx,
                                                            stim_list=stim_list,
                                                            stack_numpy=stack_numpy)
    peakresp_nmean = compute_response_peaks(dfof, peakresp_idx=peakresp_idx, stim_list=stim_list,
                                            maxk=maxk, stack_numpy=stack_numpy)
    return peakresp_nmean, baseline_avg, baseline_std


def compute_interval_response_vectors(dfof, timestamps, stim_ict, peak_start=None, peak_len=None, baseline_len=-10):
    """
    Compute response vectors, using timepoints and desired time intervals


    Parameters
    ----------
    dfof - 2D numpy array, n_cells x timestamps.size
    timestamps - 1d np.array of timestamps
    stim_ict - 1d np.array of stimulus onset times
    peak_start - length of time after stim_ict to start calculating peak_mean
    peak_len - Time period (after peak_start) over which to calculate peak_mean
               peak_mean = mean(dfof, peak_start: peak_start + peak_len)
    baseline_len - Length of time pre stim_ict over which to calculate baseline_mean and baseline_std

    Returns
    -------
    peak_mean
    baseline_mean
    baseline_std
    """
    if peak_start is None:
        peak_start = 0.5
    if peak_len is None:
        peak_len = 5
    if baseline_len is None:
        baseline_len = -10

    if dfof.shape[1] != timestamps.size:
        raise ValueError("dfof array and timestamps 1d array are not compatible sizes")

    n_trials = stim_ict.size
    n_cells = dfof.shape[0]

    peak_mean = np.empty((n_cells, n_trials))
    baseline_mean = np.empty((n_cells, n_trials))
    baseline_std = np.empty((n_cells, n_trials))

    for i, stim_on in enumerate(stim_ict):
        peak_idx = (timestamps >= (stim_on + peak_start)) \
                   & (timestamps < (stim_on + peak_start + peak_len))
        peak_mean[:, i] = np.mean(dfof[:, peak_idx], axis=1)

        baseline_idx = (timestamps >= (stim_on + baseline_len)) & (timestamps < stim_on)
        baseline_mean[:, i] = np.mean(dfof[:, baseline_idx], axis=1)
        baseline_std[:, i] = np.std(dfof[:, baseline_idx], axis=1)

    return peak_mean, baseline_mean, baseline_std



def compute_sparsity(n_std, peakresp_amp, baseline_std):
    K, n_trials = peakresp_amp.shape
    sparsity = np.empty([n_std.size, n_trials])
    for i, val in enumerate(n_std):
        thresh = val * baseline_std
        responders = peakresp_amp > (val * baseline_std)
        sparsity[i, :] = np.sum(responders, axis=0) / K
    return sparsity


def compute_sparsity_thr(thr_list, peakresp_amp, baseline_std, std_thr=3.0):
    K, n_trials = peakresp_amp.shape
    sparsity = np.empty([thr_list.size, n_trials])
    for i, val in enumerate(thr_list):
        responders = (peakresp_amp > (std_thr * baseline_std)) & (peakresp_amp > val)
        sparsity[i, :] = np.sum(responders, axis=0) / K
    return sparsity


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
