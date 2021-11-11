import numpy as np
from ryim import dirutils
import rearing_datasets, step00_process_suite2p
from rearing.ca import trace
from step00_process_suite2p import load_metadata, get_title_str, postprocess_metadata, load_all
from scipy.io import savemat


def odor_tensor_mean(cells_x_trials_x_time, df_stimulus, odor):
	"""
	For trial tensor w/ repeated trials of the same odor, calculates mean trial for odor.

	:param cells_x_trials_x_time: trial tensor w/ dimensions as indicated in name
	:type cells_x_trials_x_time: np.ndarray
	:param df_stimulus:
	:type df_stimulus: pd.DataFrame
	:param odor: which odor to find the average of
	:type odor: str
	:return:
	:rtype: np.ndarray
	"""
	odor_trials = cells_x_trials_x_time[:, df_stimulus.stimname == odor, :]
	odor_mean = odor_trials.mean(axis=1)
	return odor_mean


def align_signals(traces, timestamps, stim_ict, response_len=30, baseline_len=15):
	"""

    Takes a 2D [# cells x timepoints] array (dfof) w/ timestamps, and creates a 4D array, with axes (time, trial, cell)
	...

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
	pre_frames = int(np.ceil(baseline_len * fs))
	post_frames = int(np.ceil(response_len * fs))

	trial_T = pre_frames + post_frames
	trial_ts = (np.arange(trial_T) - pre_frames) / fs

	# empty 3D np.array
	trial_traces = np.empty([n_cells, n_trials, trial_T])
	for cid in range(n_cells):
		for i, ici in enumerate(stim_ici):
			trial_traces[cid, i, :] = traces[cid, (ici - pre_frames):(ici + post_frames)]
	return trial_traces, trial_ts


def extract(f):
	"""
    Load traces from suite2p/**/step00_process_suite2p (folder in the same directory as stat.npy)

    :param f: Path to stat.npy, uses this to parse fly information and load other dataset info
    :return: nt_meta_data, title_str, processed_traces

    """
	nt_meta_data = postprocess_metadata(load_metadata(f))
	ainfo = dirutils.parse_dataset_info(f)
	title_str = get_title_str(ainfo, nt_meta_data)

	folder = f.parent.joinpath('step00_process_suite2p')
	C_dec = np.load(folder.joinpath('C_dec.npy'))
	Fcc = np.load(folder.joinpath('Fcc.npy'))
	Fq = np.load(folder.joinpath('Fq.npy'))
	Sp = np.load(folder.joinpath('Sp.npy'))
	oasis_results = np.load(folder.joinpath('oasis_results.npy'), allow_pickle=True).item()

	processed_traces = dict(C_dec=C_dec, Fcc=Fcc, Fq=Fq, Sp=Sp, oasis_results=oasis_results)
	return nt_meta_data, title_str, processed_traces


def transform(processed_traces, timestamps, stim_ict, response_len=30, baseline_len=15):
	"""
    Takes fluorescence timeseries and converts 3D np.ndarrays [n_cells x trials x time]

    :param processed_traces: dict, should contain fluorescence timeseries and oasis_results from
                                suite2p/**/step00_process_suite2p
    :param timestamps: 1d np.ndarray
    :param stim_ict: stimulus-on times
    :param response_len: time in seconds after stimulus onset to include in trial
    :param baseline_len: time in seconds before stimulus onset to include in trial
    :return: dict, contains tensors [n_cells x trials x time] and trialtensor_ops
    """
	# Create trial tensors, trials_x_cells_x_time
	# --------------------------------------------
	tensor_C_dec, trial_ts = align_signals(processed_traces['C_dec'], timestamps, stim_ict,
	                                       response_len=response_len,
	                                       baseline_len=baseline_len)

	tensor_Sp, _ = align_signals(processed_traces['Sp'], timestamps, stim_ict,
	                             response_len=response_len,
	                             baseline_len=baseline_len)

	tensor_Fcc, _ = align_signals(processed_traces['Fcc'], timestamps, stim_ict,
	                              response_len=response_len,
	                              baseline_len=baseline_len)

	tensor_Fq, _ = align_signals(processed_traces['Fq'], timestamps, stim_ict,
	                             response_len=response_len,
	                             baseline_len=baseline_len)

	trialtensor_ops = dict(response_len=response_len, baseline_len=baseline_len)

	trial_tensors = dict(trial_ts=trial_ts,
	                     tensor_C_dec=tensor_C_dec,
	                     tensor_Sp=tensor_Sp,
	                     tensor_Fcc=tensor_Fcc,
	                     tensor_Fq=tensor_Fq,
	                     trialtensor_ops=trialtensor_ops)

	return trial_tensors


def save(outputdir, trial_tensors):
	"""
    Writes trial tensor arrays to folder.

    Should result in the following files written to outputdir:
        - trial_ts.npy
        - tensor_C_dec.npy
        - tensor_Sp.npy
        - tensor_Fcc.npy
        - trialtensor_ops.npy

    """
	# if outputdir does not exist, create it
	outputdir.mkdir(parents=True, exist_ok=True)

	# save contents of trial_tensors to .npy files
	for key, value in trial_tensors.items():
		np.save(outputdir.joinpath(f"{key}.npy"), value, allow_pickle=True)

	savemat(outputdir.joinpath("trial_tensors.mat"), trial_tensors)
	print(f'Results of step02_extract_trial_tensors saved successfully to:\n{outputdir}.')
	return True


def main(f, trialtensor_ops=None):
	if trialtensor_ops is None:
		trialtensor_ops = dict(response_len=30, baseline_len=15)

	# extract data from suite2p/**/step00_process_suite2p
	nt_meta_data, title_str, processed_traces = extract(f)

	# transform
	trial_tensors = transform(processed_traces,
	                          nt_meta_data.timestamps,
	                          nt_meta_data.df_stimulus['stim_ict'].to_numpy(),
	                          **trialtensor_ops)

	# save results to step02_extract_trial_tensors
	outputdir = f.parent.joinpath('step02_extract_trial_tensors')
	outputdir.mkdir(parents=True, exist_ok=True)
	save(outputdir, trial_tensors)
	return True


if __name__ == '__main__':
	pfo_stat_files = rearing_datasets.get_pfo_statfiles()
	odor_stat_files = rearing_datasets.get_odor_statfiles()
	stat_files = pfo_stat_files + odor_stat_files

	trialtensor_ops = dict(response_len=30, baseline_len=15)

	for f in stat_files:
		main(f, trialtensor_ops)
