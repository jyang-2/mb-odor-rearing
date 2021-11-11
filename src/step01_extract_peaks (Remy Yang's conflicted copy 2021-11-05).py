import numpy as np
from ryim import dirutils
from rearing.ca import respvec
import rearing_datasets
from step00_process_suite2p import load_metadata, get_title_str, postprocess_metadata


# %% correlation analysis, for sanity check

def compute_interval_response_vectors(dfof, timestamps, stim_ict, peak_start=None, peak_len=None, baseline_len=10):
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

        baseline_idx = (timestamps >= (stim_on - baseline_len)) & (timestamps < stim_on)
        baseline_mean[:, i] = np.mean(dfof[:, baseline_idx], axis=1)
        baseline_std[:, i] = np.std(dfof[:, baseline_idx], axis=1)
    return peak_mean, baseline_mean, baseline_std


def find_max_peaks_and_times(dfof, timestamps, stim_ict, peak_start=0.5, peak_len=5):
    """
    Finds maximum of dfof traces within time interval after stimulus onset.

    :param dfof:
    :param timestamps:
    :param stim_ict: time when stimulus is turned on
    :param peak_start:
    :param peak_len:
    :return:
    """
    max_peaks = np.zeros((dfof.shape[0], stim_ict.size))
    peak_times = np.zeros((dfof.shape[0], stim_ict.size))

    for i, stim_on in enumerate(stim_ict):
        peak_idx = (timestamps >= (stim_on + peak_start)) \
                   & (timestamps < (stim_on + peak_start + peak_len))

        dfof_peak = dfof[:, peak_idx]
        t = timestamps[peak_idx]
        max_peaks_vec = np.max(dfof_peak, axis=1)
        peak_times_vec = np.take_along_axis(t, np.expand_dims(np.argmax(dfof_peak, axis=1), axis=-1), axis=-1)

        max_peaks[:, i] = max_peaks_vec
        peak_times[:, i] = peak_times_vec

    return max_peaks, peak_times


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


def transform(dfof, timestamps, stim_ict, peak_start=0.5, peak_len=5, baseline_len=-10):
    """
    Takes fluorescence timeseries and converts to matrices with trial peak measurements, like peak_mean.

    :param dfof:
    :param timestamps:
    :param stim_ict:
    :param peak_start:
    :param peak_len:
    :param baseline_len:
    :return:
    """
    peak_mean, baseline_mean, baseline_std = compute_interval_response_vectors(dfof,
                                                                               timestamps, stim_ict,
                                                                               peak_start=peak_start,
                                                                               peak_len=peak_len,
                                                                               baseline_len=baseline_len)

    peak_amp = peak_mean - baseline_mean

    respvec_ops = dict(peak_start=peak_start, peak_len=peak_len, baseline_len=baseline_len)
    peak_data = dict(peak_mean=peak_mean, baseline_mean=baseline_mean, baseline_std=baseline_std,
                     peak_amp=peak_amp, respvec_ops=respvec_ops)
    return peak_data


def save(outputdir, peak_data):
    """
    Writes peak measurement arrays to folder.

    Should result in the following files written to outputdir:
        - peak_amp.npy
        - peak_mean.npy
        - baseline_mean.npy
        - respvec_ops.npy

    """
    outputdir.mkdir(parents=True, exist_ok=True)

    # np.save(outputdir.joinpath('peak_amp.npy'), peak_data)
    # np.save(outputdir.joinpath('peak_mean.npy'), peak_mean)
    # np.save(outputdir.joinpath('baseline_mean.npy'), baseline_mean)
    # np.save(outputdir.joinpath('respvec_ops.npy'), peak_times)

    for key, value in peak_data.items():
        np.save(outputdir.joinpath(f"{key}.npy"), value)

    print(f'Results of step01_extract_peaks saved successfully to {outputdir}.')
    return True


def main(f, respvec_ops=None):
    if respvec_ops is None:
        respvec_ops = dict(peak_start=0.5, peak_len=5, baseline_len=10)

    # extract
    nt_meta_data, title_str, C_dec = extract(f)

    stim_ict = nt_meta_data.df_stimulus['stim_ict']

    # transform
    peak_data = transform(C_dec, nt_meta_data.timestamps, stim_ict, **respvec_ops)

    # save
    outputdir = f.parent.joinpath('step01_extract_peaks')
    outputdir.mkdir(parents=True, exist_ok=True)
    save(outputdir, peak_data)
    return True


if __name__ == '__main__':
    pfo_stat_files = rearing_datasets.get_pfo_statfiles()
    odor_stat_files = rearing_datasets.get_odor_statfiles()

    respvec_ops = dict(peak_start=0.5, peak_len=5, baseline_len=-10)

    for f in pfo_stat_files:
        main(f, respvec_ops=respvec_ops)

    for f in odor_stat_files:
        main(f, respvec_ops=respvec_ops)

# %%
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


# %%

# def plot_cell(cid_):
# 	ax = draw_oasis_results(d[cid_, :], C_dec[cid_, :], Sp[cid_, :],
# 	                        sn_=Sn[cid_], x_=timestamps)
# 	ax.set_title(f"cell # {cid_} (snr={round(Sn[cid_], 2)})", fontsize=8)
# 	ax.spines[["left", "bottom"]].set_position(("data", 0))
# 	plt.subplots_adjust(bottom=0.05, top=0.9, hspace=0.25)
# 	plt.suptitle(title_str, fontsize=10)
# 	plt.show(block=False)
#
#
# # %% plot oasis results
# save_oasis_plots = True
# if save_oasis_plots:
# 	# of cells per figure
# 	n_axes = 40
# 	for cid0 in range(0, n_axes, 41):
# 		fig, axarr = plt.subplots(10, int(np.ceil(n_axes / 10)), figsize=(11, 8),
# 		                          sharey='all', sharex='all')
#
# 		for cid, ax in zip(range(cid0, cid0 + n_axes), axarr.flatten()):
# 			ax = draw_oasis_results(d[cid, :], C_dec[cid, :], Sp[cid, :],
# 			                        sn_=Sn[cid], x_=timestamps, ax_=ax)
#
# 			ax.set_title(f"cell # {cid} (snr={round(Sn[cid], 2)})", fontsize=8)
# 			ax.spines[["left", "bottom"]].set_position(("data", 0))
# 	# ax.set_ylim([-2, 10])
#
# 	plt.subplots_adjust(bottom=0.05, top=0.9, hspace=0.25)
# 	plt.suptitle(title_str, fontsize=10)
# 	plt.show(block=False)
