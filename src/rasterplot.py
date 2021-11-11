import sys
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from pathlib import Path

from scipy.stats import zscore
import math

import rastermap
from rastermap import Rastermap

from scipy.ndimage import gaussian_filter1d

suite2p_ops = {'n_components': 1, 'n_X': 100, 'alpha': 1., 'K': 1.,
               'nPC': 200, 'constraints': 2, 'annealing': True, 'init': 'pca',
               }
"""
Parameters

Rastermap first takes the specified PCs of the data, and then embeds them into n_X clusters. It returns upsampled 
cluster identities (n_X x upsamp). Clusters are also computed across Y (n_Y) and smoothed, to help with fitting. 

    n_components : int, optional (default: 2) dimension of the embedding space
    n_X : int, optional (default: 40) size of the grid on which the Fourier modes are rasterized
    nPC : nparray, int, optional (default: 400) how many of the top PCs to use during optimization
    alpha : float, optional (default: 1.0) exponent of the power law enforced on component n as: 1/(K+n)^alpha
    K : float, optional (default: 1.0) additive offset of the power law enforced on component n as: 1/(K+n)^alpha
    init : initialization of algorithm (default: 'pca') can use 'pca', 'random', or a matrix n_samples x n_components

    # hidden 'init' allowed values: 'pca', 'dbisomap', 'random', 'laplacian', 'lle'

Outputs

Rastermap model has the following attributes after running 'fit':

    embedding : array-like, shape (n_samples, n_components) Stores the embedding vectors.
    u,sv,v : singular value decomposition of data S, potentially with smoothing
    isort1 : sorting along first dimension (n_samples) of matrix
    cmap : correlation of each item with all locations in the embedding map (before upsampling)
    A : PC coefficients of each Fourier mode

"""


def suite2p_reordered_trial_display(sp, df_stimulus, frametimes, nbin=None, nplotrows=400,
                                    pre_stim=-20, post_stim=40,
                                    figure_kwargs={}, imshow_kwargs={}):

    df = df_stimulus.sort_values(by='stimname')
    n_cells = sp.shape[0]

    if nbin is None:
        # selects nbin so that the rasterplot has ~400 rows
        nbin = int(n_cells / nplotrows)
        nbin = round(nbin / 10) * 10

    img = running_average(sp, nbin)
    img = zscore(img, axis=1)

    # plot each trial
    fig_kwargs = dict(nrows=1, ncols=df_stimulus.shape[0], figsize=(24, 12))
    fig_kwargs.update(figure_kwargs)
    fig, axarr = plt.subplots(**fig_kwargs)

    im_kwargs = dict(vmin=-0.5, vmax=0.5, aspect='auto', cmap='gray_r', origin='lower')
    im_kwargs.update(imshow_kwargs)

    for odor, t_stim, ax in zip(df['stimname'], df['stim_ict'], axarr.flat):
        # time intervals, before and after stimulus onset
        # t = [start, stim onset, end] for trial
        t = [t_stim + pre_stim, t_stim, t_stim + post_stim]

        idx = np.floor(np.interp(t, frametimes, np.arange(frametimes.size))).astype(int)
        print(idx, type(idx[0]))

        im = ax.imshow(img[:, idx[0]:idx[-1]], **im_kwargs)
        ax.set_title(f"{odor}")
        ax.axvline(idx[1] - idx[0], linestyle='--', linewidth=1)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    axarr[0].get_yaxis().set_visible(True)
    fig.colorbar(im)

    return fig, axarr


    fig.suptitle(f"{title_str}\nconc = {blk_conc}")



def norm_traces(traces):
    traces = np.squeeze(traces)
    traces = zscore(traces, axis=1)
    traces = np.maximum(-4, np.minimum(8, traces)) + 4
    traces /= 12
    return traces


def prep_to_plot(spF):
    spF = zscore(spF, axis=1)
    spF = np.minimum(8, spF)
    spF = np.maximum(-4, spF) + 4
    spF /= 12
    return spF


# this function performs a running average filter over the first dimension of X
# (faster than doing gaussian filtering)
def running_average(x, nbin=50):
    Y = np.cumsum(x, axis=0)
    Y = Y[nbin:, :] - Y[:-nbin, :]
    return Y


def suite2p_display(S, isort=None, nbin=50, figure_kwargs={}, imshow_kwargs={}):
    if isort is not None:
        S = S[isort, :]

    Sfilt = running_average(S, nbin)
    Sfilt = zscore(Sfilt, axis=1)

    figure_kwargs0 = dict(figsize=(12, 6))
    figure_kwargs0.update(figure_kwargs)

    imshow_kwargs0 = dict(vmin=-0.5, vmax=5, aspect='auto', cmap='gray_r', origin='lower')
    imshow_kwargs0.update(imshow_kwargs)

    fig1, ax = plt.subplots(**figure_kwargs0)
    img = ax.imshow(Sfilt, **imshow_kwargs0)
    fig1.colorbar(img)

    plt.xlabel('time points')
    plt.ylabel('sorted neurons')
    return fig1, ax


def main(statfile, save_outputs=True, nbin=None, figure_kwargs={}, imshow_kwargs={}):
    if isinstance(statfile, str):
        statfile = Path(statfile)

    stat = np.load(statfile, allow_pickle=True)[0]
    iscell = np.load(statfile.with_name('iscell.npy'))
    iscell = iscell[:, 0].astype(np.bool_)
    Fcell = np.load(statfile.with_name('F.npy'))
    Fneu = np.load(statfile.with_name('Fneu.npy'))
    F = Fcell - Fneu * 0.7

    # select good cells
    sp = F[iscell, :]
    spn = norm_traces(sp)

    if nbin is None:
        nbin = math.ceil(sp.shape[0] / 400)
        print(f"nbin = {nbin}")

    use_default = True
    if use_default:
        ops = suite2p_ops
    else:
        ops = {'n_components': 2,  # (default: 2) dimension of the embedding space
               'n_X': 100,  # (default: 40) size of the grid on which the Fourier modes are rasterized
               'alpha': 1.0,  # (default: 1.0) exponent of the power law enforced on component n as: 1/(K+n)^alpha
               'K': 1.0,  # (default: 1.0) additive offset of the power law enforced on component n as: 1/(K+n)^alpha
               'nPC': 200,  # (default: 400) how many of the top PCs to use during optimization
               'constraints': 2,  # (default: 1.0) exponent of the power law enforced on component n as: 1/(K+n)^alpha
               'annealing': True,
               'init': 'pca',  # (default: 'pca') can use 'pca', 'random', or a matrix n_samples x n_components
               }
    model = Rastermap(**ops)

    # fit does not return anything, it adds attributes to model
    # attributes: embedding, u, s, v, isort1
    embedding = model.fit_transform(spn)
    isort1 = np.argsort(embedding[:, 0])

    # %% plot rastermap
    spp = prep_to_plot(spn)
    fig, img = suite2p_display(spp, isort=isort1, nbin=nbin, figure_kwargs=figure_kwargs, imshow_kwargs=imshow_kwargs)

    title_str = str(Path(statfile).relative_to(Path.cwd().parent))

    plt.title(f"{title_str}\n# cells={sp.shape[0]}")
    plt.show()

    if save_outputs:
        np.save(statfile.with_name('embedding.npy'), model)
        fig.savefig(statfile.with_name('rastermap.png'))
        fig.savefig(statfile.with_name('rastermap.pdf'))

    return model, fig


# %%
if __name__ == 'main':
    folder = Path.cwd().joinpath('data', 'processed_data', '2021-08-24') \
        .rglob('**/suite2p/combined/stat.npy')
    folder = sorted(list(folder))
    for file in folder:
        main(file)
# mdl, fig_rastermap = main(sys.argv[1])
