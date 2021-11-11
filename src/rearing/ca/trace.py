import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import filters
from scipy.stats import zscore, mode
from scipy import stats, signal


# from suite2p.extraction import dcnv


def split_traces_to_blocks(F_traces, block_idx):
	"""
	splits traces into blocks (2D numpy array to a list of 2D numpy arrays)

	:param traces: [# cells x # timestamps]  Array of fluorescence traces
	:type traces: np.array
	:param block_idx: gives block index for each frame in traces
	:type: list or np.array
	:return: list of np.array fluorescence traces (split so that each item is one acquisition)

	Parameters
	----------
	F_traces
	F_traces
	"""
	blocks = np.unique(block_idx)

	blocked_traces = map(lambda x: F_traces[:, block_idx == x], blocks)
	return list(blocked_traces)


def extract_dfof(traces, fs=6.1, win=45, percentile=20):
	baseline = percentile_baseline(traces, fs=fs, win=win, percentile=percentile)
	df = traces - baseline
	dfof = df / baseline
	return dfof, df, baseline


# def suite2p_extract_dfof(F, fs, **kwargs):
# 	"""
#
#     Parameters
#     ----------
#     F : Neuropil-corrected traces, [# cells * timestamps]
#     baseline : baseline fluorescence, calculated with percentile_baseline(...)
#     win_baseline :
#     sig_baseline :
#     fs :
#     prctile_baseline :
#
#     Returns
#     -------
#
#     """
# 	F_corrected = dcnv.preprocess(
# 		F=F,
# 		fs=fs,
# 		baseline=kwargs['baseline'],
# 		win_baseline=kwargs['win_baseline'],
# 		sig_baseline=kwargs['sig_baseline'],
# 		prctile_baseline=kwargs['prctile_baseline']
# 	)
#
# 	baseline = F - F_corrected
# 	dfof = df / baseline
# 	return dfof, df, baseline


def percentile_baseline(traces, fs=1, win=60, percentile=20):
	"""
	Function to find baseline for fluorescence traces w/ running percentile window

	Example -
	-----------------
	`
	baseline = percentile_baseline(Fc, fs=6.1, win=45, percentile=10)
	`
	-----------------

	:param traces: [neurons x time] Fluorescence activity matrix, 1 column per timepoint
	:type traces: numpy.ndarray
	:param fs: Sampling rate (Hz)
	:type fs: float
	:param win: Window length of percentile filter (in seconds)
	:type win: int or float
	:param percentile: Percentile of activity to use as baseline, between 0 and <100
	:type percentile: int
	:return:
		baseline - [neurons x time] Baseline fluorescence
	:rtype: numpy.ndarray
	"""
	win_frames = int(win * fs)
	Flow = filters.percentile_filter(traces, percentile=percentile, size=(1, win_frames))
	return Flow


def preprocess_for_rastermap(traces, method='suite2p'):
	"""
	Function to zscore and normalize cell responses (from suite2p rastermap, for visualization
	purposes)
	---------------
	Example -
	`
	dfof_z = preprocess_for_rastermap(dfof)
	`
	---------------
	:param traces: [neurons x time] Fluorescence activity matrix, 1 column per timepoint
	:type traces: numpy.ndarray
	:param method: {'suite2p, }
	:type method: str
	:return:
		traces - [neurons x time], in range (0, 1)
	:rtype: numpy.ndarray
	"""
	if method == 'suite2p':
		traces = np.squeeze(traces)
		traces = zscore(traces, axis=1)
		traces = np.maximum(-4, np.minimum(8, traces)) + 4
		traces /= 12
		return traces


# design lowpass filter
def get_lowpass_filter(filt_ops):
	"""

	Parameters
	----------
	filt_ops :

	Returns
	-------

	"""
	fs = filt_ops['fs']
	cutoff = filt_ops['cutoff']
	trans_width = filt_ops['trans_width']
	numtaps = filt_ops['numtaps']
	b = signal.remez(numtaps, [0, cutoff, cutoff + trans_width, 0.5 * fs], [1, 0], Hz=fs)
	return b


def plot_filter_response(fs, w, h, title):
	"""
	Utility function to plot response functions

	example -
	---------------
	`
	fs = 6.1  # Sample rate, Hz
	cutoff = .5  # Desired cutoff frequency, Hz
	trans_width = .1  # Width of transition from pass band to stop band, Hz
	numtaps = 203  # Size of the FIR filter.
	b = signal.remez(numtaps, [0, cutoff, cutoff + trans_width, 0.5 * fs], [1, 0], Hz=fs)
	w, h = signal.freqz(b, [1], 2000)
	plot_filter_response(fs, w, h, "Low-pass Filter")
	`
	:param fs: sampling rate
	:type fs: float
	:param w: The frequencies at which h was computed, in the same units as fs.
			  By default, w is normalized to the range [0, pi) (radians/sample).
	:type w: numpy.ndarray
	:param h: The frequency response, as complex numbers.
	:type h: numpy.ndarray
	:param title: Plot title
	:type title: str
	:return:
			 **fig** - matplotlib figure handle
			  **ax** - matplotlib axes handle
	:rtype:
			matplotlib.figure.Figure
			matplotlib.axes.Axes
	"""
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(0.5 * fs * w / np.pi, 20 * np.log10(np.abs(h)))
	# ax.set_ylim(-40, 5)
	# ax.set_xlim(0, 0.5 * fs)
	ax.grid(True)
	ax.set_xlabel('Frequency (Hz)')
	ax.set_ylabel('Gain (dB)')
	ax.set_title(title)
	plt.show()
	return fig, ax


# compute framerate from timestamps


def compute_fs(frame_times, method='accurate'):
	"""
	Given a vector of timestamps (frame times), computes sampling rate fs.

	:param frame_times: List of vectors ( fti.frame_times )
	:type frame_times: np.array or List
	:param method: {'round', 'accurate'}; use 'accurate' if you plan on using nitime
	:type method: str
	:return: fs
	:rtype: float
	"""

	if method == 'round':
		frametimes = np.round(frame_times, 4)
		dt = mode(np.diff(frametimes)).mode[0]
		fs = float(round(1.0/dt, 1))
	else:
		# remove outliers/inter-acquisition dt values
		dt = np.diff(frame_times)
		dt = stats.trim1(dt, 0.01, 'right').mean()
		fs = 1/dt
	return fs
