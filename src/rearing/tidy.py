import numpy as np
import pandas as pd


def tidy_traces(dfof, df_frametimeinfo, cell_ids=None):
    """
    Tidies dfof numpy array, returns DataFrame with same columns as df_frametimeinfo, plus "cid" and "dfof"

    Used for plotting stacked calcium traces.

    Example of tidy_dfof:

    |    |   frame_num |   frame_times |   scope_pulse_idx |   olf_pulse_idx |   baseline_idx |   peakresp_idx |   trialtrace_idx |   cid |      dfof |
    |---:|------------:|--------------:|------------------:|----------------:|---------------:|---------------:|-----------------:|------:|----------:|
    |  0 |           1 |       16.2518 |                 1 |               0 |              0 |              0 |                0 |     0 | 0.319292  |
    |  1 |           2 |       16.4497 |                 1 |               0 |              0 |              0 |                0 |     0 | 0.234497  |
    |  2 |           3 |       16.6476 |                 1 |               0 |              0 |              0 |                0 |     0 | 0.163523  |
    |  3 |           4 |       16.8455 |                 1 |               0 |              0 |              0 |                0 |     0 | 0.116671  |
    |  4 |           5 |       17.0434 |                 1 |               0 |              0 |              0 |                0 |     0 | 0.0965842 |


    Parameters
    ----------
    dfof : np.ndarray, 2D
        array of neuron activity traces, w/ row=neuron and column=timepoint

    df_frametimeinfo : pd.DataFrame

    cell_ids : list, int
        If not given, cids assigned starting from 0-index

    Returns
    -------
    tidy_dfof :  pd.DataFrame
        "tidy" activity traces (melted dfof

    """
    n_cells, n_ts = dfof.shape

    if cell_ids is None:
        cids = list(range(n_cells))
    else:
        cids = cell_ids

    df_dfof = pd.DataFrame(dfof.T, columns=cids).join(df_frametimeinfo)

    tidy_dfof = pd.melt(df_dfof,
                        id_vars=df_frametimeinfo.columns,
                        value_vars=cids,
                        var_name='cid',
                        value_name='dfof'
                        )
    return tidy_dfof


def tidy_response_peaks(peakresp_max, peakmax_times, bin_response=None, df_frametimeinfo=None):
    n_cells, n_trials = peakresp_max.shape
    cid_grid, trial_grid = np.mgrid[0:n_cells, 0:n_trials]

    if bin_response is None:
        bin_response = np.ones(peakresp_max.shape, dtype=bool)

    tidy_resppeaks = pd.DataFrame(
        np.stack([cid_grid[bin_response], peakresp_max[bin_response], peakmax_times[bin_response]], axis=-1),
        columns=['cid', 'peakresp_max', 'peakmax_time'])
    tidy_resppeaks['cid'] = tidy_resppeaks['cid'].astype(int)

    if df_frametimeinfo is not None:
        tidy_resppeaks = tidy_resppeaks.join(
            df_frametimeinfo.set_index('frame_times').loc[tidy_resppeaks['peakmax_time']].reset_index())

    return tidy_resppeaks


def tidy_df_corrmat(df):
    tidy_df = df.melt(ignore_index=False) \
        .reset_index() \
        .dropna()
    tidy_df['stimname_col'] = tidy_df['stimname_col'].astype(tidy_df['stimname_row'].dtype)
    return tidy_df


def tidy_triu_df_corrmat(df):
    mask = np.triu(np.ones_like(df, dtype=np.bool))
    tidy_triu = df.mask(mask) \
        .melt(ignore_index=False) \
        .reset_index() \
        .dropna()
    tidy_triu['stimname_col'] = tidy_triu['stimname_col'].astype(tidy_triu['stimname_row'].dtype)
    # grouped = tidy_triu.groupby(['stimname_row', 'stimname_col'], observed=True, as_index=False)

    return tidy_triu


def tidy_tril_df_corrmat(df):
    mask = np.tril(np.ones_like(df, dtype=np.bool))
    tidy_tril = df.mask(mask) \
        .melt(ignore_index=False) \
        .reset_index() \
        .dropna()
    tidy_tril['stimname_col'] = tidy_tril['stimname_col'].astype(tidy_tril['stimname_row'].dtype)
    # grouped = tidy_triu.groupby(['stimname_row', 'stimname_col'], observed=True, as_index=False)

    return tidy_tril

#def tidy_dfof
