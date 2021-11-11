import numpy as np
import pandas as pd
import xmltodict


def compute_peakresp_idx(frame_times, stim_ict, peakresp_win):
    """
    Computes peakresp_idx (column in df_frametimeinfo).

    Has same length as frame_times, but nonzero indices indicate time windows for response peak calculations.


    Parameters
    ----------
    frame_times : float, 1D numpy array
        timestamps for every movie frame

    stim_ict : float, list
        stimulus start time (in seconds)

    peakresp_win : float
        window (in seconds) after stim_ict for computing response amplitudes


    Returns
    -------
    peakresp_idx : float, 1D array
        size [timepoints,] >0 values for indexing into response windows

        example: [0, 0, ... 0, 1, 1, 1, 1, 0, 0, 0, ... 0, 2, 2, 2, 2, 0, 0, 0]
    """

    n_stim = len(stim_ict)
    peakresp_idx = np.zeros(frame_times.shape, dtype=int)

    stim_list = np.arange(n_stim) + 1

    for i, stim_on in enumerate(stim_ict):
        peakresp_idx[(stim_on <= frame_times) & (frame_times < (stim_on + peakresp_win))] = i + 1

    return peakresp_idx
