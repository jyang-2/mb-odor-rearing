import numpy as np
from scipy.signal import find_peaks

from rearing.ca import trace

"""
2D color video:  (t, row, col, ch)

3D color video: (t, z, row, col, ch) 
"""


def align_signals(traces, timestamps, stim_ict, response_len=30, baseline_len=15):
	"""

	Takes a 2D [# cells x timepoints] array (dfof) w/ timestamps, and creates a 4D array, with axes (time, trial, cell)
	(cell, trial, time)

	Parameters
	----------
	response_len
	traces
	timestamps
	stim_ict
	baseline_len

	Returns
	-------
		trial_traces = 3D np.array, [cells x trials x timepoints]

		trial_ts = 1D np.array, timestamps for trial timepoints, w/ t=0 indicating stimulus onset.
	"""

	n_cells, T = traces.shape
	n_trials = stim_ict.size

	fs = trace.compute_fs(timestamps)

	baseline_ict = stim_ict - baseline_len  # initial time
	response_fct = stim_ict + response_len  # final (end of response trace) time

	stim_ici = np.rint(np.interp(stim_ict, timestamps, np.arange(timestamps.size))).astype(int)
	pre_frames = int(np.ceil(baseline_len*fs))
	post_frames = int(np.ceil(response_len*fs))

	trial_T = pre_frames + post_frames
	trial_ts = (np.arange(trial_T) - pre_frames) / fs

	# empty 3D np.array
	trial_traces = np.empty([n_cells, n_trials, trial_T])
	for cid in range(n_cells):
		for i, ici in enumerate(stim_ici):
			trial_traces[cid, i, :] = traces[cid, (ici-pre_frames):(ici+post_frames)]
	return trial_traces, trial_ts



