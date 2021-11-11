import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import seaborn as sns
from pathlib import Path
import yaml
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.subplots as sp
import re
import os
import pprint
# from drosolf import orns, pns
from scipy import signal
import scipy.io as sio
from scipy.stats import zscore
import importlib
import textwrap
from itertools import product
from rearing import hallem, io, meta
from scipy.ndimage import filters
from matplotlib.pyplot import gcf
from mpl_toolkits.axes_grid1 import make_axes_locatable
from rearing.ca import trace
import itertools

import rastermap

plt.rcParams.update({'font.size': 10, 'backend': 'Qt5Agg'})
# plt.rcParams['backend'] = 'module://backend_interagg'


def blocked_stackplot(F_traces, timestamps, cids=None, df_stimulus=None, df_frametimeinfo=None,
                      peak_start=None, peak_len=None, stim_to_plot=None, F_overlay=None):
    """

    Parameters
    ----------
    F_traces
    timestamps
    cids
    df_stimulus
    df_frametimeinfo
    peak_start
    peak_len
    stim_to_plot

    Returns
    -------
        fig, axarr
    """
    # plot first 20 cells if cids is not given
    if cids is None:
        if isinstance(F_traces, list):
            cids = list(range(np.min(F_traces[0].shape[0], 20)))
        else:
            cids = list(range(np.min([F_traces.shape[0], 20])))
    # of acquisitions

    if isinstance(F_traces, np.ndarray):    # split ndarray using df_frametimeinfo
        F_traces = trace.split_traces_to_blocks(F_traces, df_frametimeinfo.scope_pulse_idx)

    if F_overlay is not None:
        if isinstance(F_overlay, np.ndarray):
            F_overlay = trace.split_traces_to_blocks(F_overlay, df_frametimeinfo.scope_pulse_idx)
        draw_overlay = True
    else:
        draw_overlay = False

    nacq = len(F_traces)
    n_cells = F_traces[0].shape[0]
    scope_pulse_idx = df_frametimeinfo['scope_pulse_idx'] - 1

    color_dict = {'1-6ol': 'tab:blue',
                  '1-5ol': 'tab:green',
                  'EP': 'tab:red',
                  '1-6ol+EP': 'tab:purple',
                  'VA': 'tab:brown',
                  'pfo': 'tab:gray'}

    if stim_to_plot is None:
        stim_to_plot = color_dict.keys()

    uspi = scope_pulse_idx.unique().tolist()

    # plot fluorescence
    fig, axarr = plt.subplots(len(cids), nacq, figsize=(12, 12), sharey=True, sharex='col')
    for cid_row, cid in enumerate(cids):
        if cid < n_cells:
            axarr[cid_row, 0].set(ylabel=str(cid))
            for i in scope_pulse_idx.unique():
                ts = timestamps[scope_pulse_idx == i]

                if draw_overlay:
                    axarr[cid_row, uspi.index(i)].plot(ts, F_traces[uspi.index(i)][cid, :], label='dfof', color='k',
                                                       linewidth=0.5, alpha=0.5)
                    axarr[cid_row, uspi.index(i)].plot(ts, F_overlay[uspi.index(i)][cid, :], label='deconv.', color='k',
                                                       linewidth=0.5)
                else:
                    axarr[cid_row, uspi.index(i)].plot(ts, F_traces[uspi.index(i)][cid, :], label='dfof', color='k',
                                                       linewidth=0.5)
    for ax in axarr.flat:
        ax.label_outer()

    # add stimulus lines
    acq_idx = np.interp(df_stimulus['stim_ict'],
                        df_frametimeinfo['frame_times'],
                        df_frametimeinfo['scope_pulse_idx']).astype('int')

    _, _, acq_idx = np.unique(acq_idx, return_index=True, return_inverse=True)

    for idx, col, ostim in zip(df_stimulus['stim_ict'], acq_idx, df_stimulus['stimname']):  #
        for cid_row, cid in enumerate(cids):
            if ostim in stim_to_plot:

                if (peak_start is not None) and (peak_len is not None):  # plot vspan
                    axarr[cid_row, col].axvspan(idx + peak_start, idx + peak_start + peak_len,
                                            label=ostim,
                                            alpha=0.3, color=color_dict[ostim])
                else:  # plot line
                    axarr[cid_row, col].axvline(idx, label=ostim, ymin=-10, ymax=1,
                                            linewidth=1, linestyle='-', alpha=0.3, color=color_dict[ostim])

    plt.legend(bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, fontsize=8)

    return fig, axarr


def plot_rastermap(sp, isort, title='Neural data embedding', clim=(0.3, 0.7), box_aspect=1,
                   stimframes=None, blockframes=None, xid=None):
    splot = sp[isort, :]

    # fig, (ax0, ax) = plt.subplots(ncols=2, nrows=1, sharey='row', figsize=(11, 4))
    #
    # gs_kw = dict(height_ratios=[1], width_ratios=[1, 100])
    # fig, (ax0, ax) = plt.subplots(ncols=2, nrows=1, constrained_layout=True,
    #                               gridspec_kw=gs_kw, sharey='row', figsize=(11, 6))

    fig = plt.figure(figsize=(11, 6))
    axim = plt.imshow(splot, vmin=clim[0], vmax=clim[1], cmap='gray_r', interpolation='none', aspect='auto')
    # axim = ax.imshow(splot, vmin=clim[0], vmax=clim[1], cmap='gray_r', interpolation='none', aspect='auto')
    ax = plt.gca()
    # ax.set_box_aspect(box_aspect)

    linealpha = 1.0
    linewidth = 0.5
    if stimframes is not None:
        [plt.axvline(x=x, linewidth=linewidth, alpha=linealpha) for x in stimframes]

    if blockframes is not None:
        [plt.axvline(x=x, color='k', linewidth=linewidth, alpha=linealpha) for x in blockframes]

    plt.xticks(rotation=45, fontsize=8, ha='right')
    plt.yticks(fontsize=8)
    plt.ylabel('neurons')
    plt.xlabel('time points (frames)')
    plt.suptitle(title)
    fig.colorbar(axim, ax=ax, shrink=0.4)

    if xid is not None:
        cclust = xid[isort]
        lloc = np.argwhere(np.diff(cclust)) + 0.5
        [plt.axhline(y=y, color='k', linewidth=linewidth, alpha=linealpha) for y in lloc]

    # cclust = np.array([item % 2 for item in cclust])
    # cclust = cclust[:, np.newaxis]
    # ax0.imshow(cclust, cmap='gray_r', vmin=0, vmax=1, aspect='equal')

    # fig.set_constrained_layout_pads(w_pad=0, h_pad=0, hspace=0, wspace=0)
    plt.show()

    return fig, axim


def ortho_pcolormesh(activity, scope_pulse_idx=None, olf_pulse_idx=None, cell_clusters=None):
    K, T = activity.shape

    fig, ax = plt.subplots(figsize=(11, 6))

    y = np.arange(1, K + 1)
    x = np.arange(T)

    divider = make_axes_locatable(ax)
    # pos = ax.get_position()

    ax_x_scope = divider.append_axes("top", size="4%", pad="1%")
    ax_x_olf = divider.append_axes("top", size="4%", pad="1%")
    ax_y = divider.append_axes("right", size="2%", pad="1%")

    ax.set_ylabel('cells')
    ax.set_xlabel('frames')

    # scope pulse idx (top)
    ax_x_scope.xaxis.set_visible(False)
    ax_x_scope.set_yticks([0.5])
    ax_x_scope.set_yticklabels(['scope_pulse_idx'], rotation='horizontal',
                               horizontalalignment='right', verticalalignment='center', fontsize=8)
    # olf pulse idx (top)
    ax_x_olf.xaxis.set_visible(False)
    ax_x_olf.set_yticks([0.5])
    ax_x_olf.set_yticklabels(['olf_pulse_idx'], rotation='horizontal',
                             horizontalalignment='right', verticalalignment='center', fontsize=8)
    # cell cluster idx (right)
    ax_y.set_xticks([])
    ax_y.set_yticks([])
    ax_y.yaxis.set_label_position('right')
    ax_y.set_ylabel('cell clusters', color='black', fontsize=8)

    cf = ax.pcolormesh(x, y, activity, shading='auto', vmin=0.0, vmax=1.0, cmap='rocket', rasterized=True)
    if scope_pulse_idx is not None:
        ax_x_scope.pcolormesh(scope_pulse_idx, cmap='tab20c')
    if olf_pulse_idx is not None:
        ax_x_olf.pcolormesh(olf_pulse_idx, cmap='tab20c')
    if cell_clusters is not None:
        ax_y.pcolormesh(cell_clusters, cmap='prism')

    return fig, [ax, [ax_x_scope, ax_x_olf], ax_y]


def make_sparsity_plot(df, stim_order=None, title=None, threshold_type='n_std'):
    """
   Makes facetgrid of sparsities plotted as a function of n_std threshold for peak amplitude response.


    :param threshold_type:
    :type threshold_type:
    :param df: SparsityAnalysis.df_melted
    :type df: pd.DataFrame
    :param stim_order: Ordering of desired stimulus types (row)
    :param title:  suptitle for plot
    :type title: str
    :type stim_order: list
    :return: g
    :rtype: sns.FacetGrid
    """

    sns.set_style('darkgrid')
    g = sns.FacetGrid(df, col="stimname", hue='conc', col_order=stim_order,
                      aspect=0.8, height=4, sharex=True, sharey=True)

    g.map(sns.lineplot, "n_std", "sparsity")
    g.add_legend()
    g.map(plt.axhline, y=0.1, ls="-", c=".5", alpha=0.5, zorder=1)
    g.map(plt.axhline, y=0.05, ls="-", c=".5", alpha=0.5, zorder=1)
    g.set_axis_labels(threshold_type, "response sparsity")

    g.fig.subplots_adjust(top=0.8, left=0.1, bottom=0.1)  # adjust the Figure in rp
    if title is not None:
        g.fig.suptitle(title, fontsize=12)
    return g


def add_trace(ax, x, y, cell_peakmax_time=None, cell_peakmax=None):
    lh = ax.plot(x, y, linewidth=1)

    if (cell_peakmax_time is not None) and (cell_peakmax is not None):
        pm_idx = np.nonzero((cell_peakmax_time >= x[0]) and (cell_peakmax_time < x[-1]))
        ax.plot(cell_peakmax_time[pm_idx], cell_peakmax[pm_idx], linestyle=None, marker=11)
    return lh


def add_peak_markers(ax, cell_peakmax_time, cell_peakmax, markershift=0.1):
    mh = ax.plot(cell_peakmax_time, cell_peakmax + markershift, linestyle="None", marker=11, color='r')
    return mh


def transform_trace(y, cid, offset=1.0, scaling=1.0):
    ys = y * scaling
    cid_offset = cid * offset
    yt = ys - cid_offset
    return yt


def stacked_traces(traces, offset=1.0, scaling=1.0, df_stimulus=None, df_frametimeinfo=None,
                   peakmax_idx=None, bin_response=None):
    """

    :param traces:
    :param offset:
    :param scaling:
    :param df_stimulus:
    :param df_frametimeinfo:
    :param peakmax_idx:
    :param bin_response:
    :return:
    """
    n_cells, n_timepoints = traces.shape

    if df_frametimeinfo is not None:
        block_traces = trace.split_traces_to_blocks(traces, df_frametimeinfo['scope_pulse_idx'])
        n_blocks = len(block_traces)
    else:
        block_traces = [traces]
        n_blocks = 1

    fig, axarr = plt.subplots(1, n_blocks, figsize=(7, 14), sharey=True)

    cell_locs = [None] * n_cells
    for i in range(0, n_blocks):
        y = block_traces[i]

        df = df_frametimeinfo.loc[df_frametimeinfo['scope_pulse_idx'] == i + 1]
        ts = df['frame_times'].to_numpy()

        block_start = ts[0]
        block_end = ts[-1]
        # df['peakresp_idx'].isin()

        for cid in range(n_cells):
            yt = transform_trace(y, cid, offset=offset, scaling=scaling)

            axarr[i].plot(ts, transform_trace(y, cid), linewidth=1)
            cell_locs[cid] = -cid * offset

        ylabels = [str(item) for item in range(n_cells)]
        print(len(ylabels) == len(cell_locs))
        axarr[i].set_yticks(cell_locs, minor=False)
        axarr[i].set_yticklabels(ylabels, fontsize=8)
        # draw shaded stimulus
        peakresp_idx = df['peakresp_idx'].to_numpy()
        stim_idx = np.unique(peakresp_idx[peakresp_idx > 0])

        for stim in stim_idx:
            stim_on = ts.iloc[np.argwhere(peakresp_idx == stim)[0][0]]
            stim_off = ts.iloc[np.argwhere(peakresp_idx == stim)[-1][0]]
            axarr[i].axvspan(stim_on, stim_off, facecolor='k', alpha=0.3)
            axarr[i].text()

    return fig, axarr


def mean_corrmats(dist, annot=True):
    """
    Makes figure w/ various mean correlation heatmaps (cosine, pearson, spearman, kendall).

    Parameters
    ----------
    dist :
    annot :

    Returns
    -------

    """
    fig, axarr = plt.subplots(2, 2, figsize=(8, 8), constrained_layout=True)
    # fig, axarr = plt.subplots(2, 2, tight_layout=True)
    haxarr = [None] * 4

    for i, corrtype in enumerate(['cosine', 'pearson', 'spearman', 'kendall']):
        haxarr[i] = sns.heatmap(dist[corrtype], cmap='RdBu_r', vmin=-1, vmax=1,
                                ax=axarr.flat[i],
                                annot=annot, annot_kws=dict(fontsize=8),
                                square=True,
                                linewidth=0.25,
                                cbar_kws={'shrink': .3})
        axarr.flat[i].set_xlabel('')
        axarr.flat[i].set_ylabel('')

        plt.setp(axarr.flat[i].get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor", fontsize=9)
        plt.setp(axarr.flat[i].get_yticklabels(), rotation=0, ha="right",
                 rotation_mode="anchor", fontsize=9)
        axarr.flat[i].set_title(corrtype, fontsize=10)
    return fig, axarr, haxarr


def corrmats(dist, ticklabels=None, include_trialnum=True, annot=False):
    """
    Makes figure w/ various correlation heatmaps (cosine, pearson, spearman, kendall).

     example :
    --------

    :param dist: contains dataframes w/ population response correlations
    :type dist: dict

    :param ticklabels:
    :type ticklabels:  list

    :param include_trialnum: whether or not to include trial # in ticklabels
    :type include_trialnum: bool

    :return: fig1, axarr, haxarr
    :rtype:

    Parameters
    ----------
    annot
    """
    df_corr = dist['pearson']
    if ticklabels is None:
        if include_trialnum:
            ticklabels = [f"{item[1]} ({item[0]})" for item in df_corr.columns.to_list()]
        else:
            # ticklabels = [item[1] for item in df_corr.columns.to_list()]
            ticklabels = [f"{row.Index[1]} [{row.Index[2]}]" for row in df_corr.itertuples()]

    clim = [-1.0, 1.0]
    fig1, axarr = plt.subplots(2, 2, figsize=(8, 8), constrained_layout=True)
    axarr = axarr.flatten()
    haxarr = [None] * 4

    for i, m in enumerate(['cosine', 'pearson', 'spearman', 'kendall']):
        haxarr[i] = sns.heatmap(dist[m].values, cmap='RdBu_r', vmin=clim[0], vmax=clim[1],
                                ax=axarr[i],
                                xticklabels=ticklabels, yticklabels=ticklabels,
                                annot=annot,
                                square=True, linewidths=.25, cbar_kws={"shrink": .3})

        axarr[i].set_title(m)
    return fig1, axarr, haxarr


def corrscatter(df, monotonic=False):
    """
    Takes upper triangular elements from a corrmat and plots trial-trial correlations
    (seaborn, stripplot).

    Parameters
    ----------
    df : correlation  matrix dataframe
    monotonic : whether to sort odors by ascending value

    Returns
    -------
        fig : figure handle to seaborn plot
        ax : axis handle to scatter plot of odor correlation values (upper triangular) in df

    """
    # make correlation dataframe wide to long
    keep = np.triu(np.ones(df.shape), k=1).astype('bool').reshape(df.size)

    if ('stim_idx' in df.index.names) & ('stim_idx' in df.columns.names):
        temp = df.droplevel('stim_idx', axis=1).droplevel('stim_idx', axis=0)
    elif ('stim_idx_row' in df.index.names) & ('stim_idx_col' in df.columns.names):
        temp = df.droplevel('stim_idx_col', axis=1).droplevel('stim_idx_row', axis=0)

    # dfs should be a long dataframe w/ indices ['stimname_row', 'stimname_col'], and a column named 'corr'
    dfs = temp.stack()[keep].rename('corr').to_frame()

    dfs.index.set_names(names=['stimname_a', 'stimname_b'], inplace=True)
    dfs.reset_index(inplace=True)
    dfs['stimname'] = list(zip(dfs['stimname_a'], dfs['stimname_b']))
    dfs_meansorted = dfs.loc[:, ['stimname', 'corr']].groupby(['stimname']).mean().sort_values(by='corr')

    if monotonic:
        order = dfs_meansorted.index.to_list()
    else:
        order = dfs['stimname'].unique().tolist()

    fig, ax = plt.subplots(figsize=(4, 5), constrained_layout=True)
    sns.stripplot(data=dfs, x="stimname", y="corr", ax=ax, dodge=True, order=order,
                  edgecolor='white')

    ax.set_ylim((-0.2, 1.0))
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8, rotation=45, ha="right",
                       rotation_mode="anchor")
    sns.despine()
    return fig, ax


def draw_stimuli(ax, stim_on, stim_off, labels=None):
    """
    Draw shaded background for stimulus times and label w/ stimulus names.

    Parameters
    ----------
    ax : matplotlib axes

    stim_on : float, list
        stim on time

    stim_off : float, list
        stim off time

    labels : str, list
        stimulus labels

    Returns
    -------
    ax : matplotlib axes
        axes w/ shaded stimulus times & stimulus names

        example:
            stim_on = [20, 65]
            stim_off = [22, 67]
            stimnames = ['odor 1', 'odor 2']
            ax = vis.draw_stimuli(ax, stim_on, stim_off, stimnames)
    """
    ymin, ymax = ax.get_ylim()
    dy = ymax - ymin
    ydisplay = 0.9 * dy + ymin

    if labels is None:
        labels = [''] * len(stim_on)

    for stim in zip(stim_on, stim_off, labels):
        ax.axvspan(stim[0], stim[1], facecolor='0.6')
        ax.text(stim[0], ydisplay, stim[2], ha="left", va="center", size=8)
    return ax


def plotly_correlograms(dist, title=None):
    fig = sp.make_subplots(rows=2,
                           cols=2,
                           subplot_titles=list(dist.keys()))
    sidx = {'pearson': (1, 1),
            'spearman': (1, 2),
            'kendall': (2, 1),
            'cosine': (2, 2)}
    for key, df in dist.items():
        print(key)
        fig.add_trace(
            go.Heatmap(
                z=df.values,
                x=list(map(lambda x: f"{x[1]} {x[0]:>3}", df.columns.to_list())),
                y=list(map(lambda x: f"{x[1]} {x[0]:>3}", df.index.to_list())),
                colorscale='rdbu_r',
                zmin=-1.0,
                zmax=1.0
            ),
            row=sidx[key][0],
            col=sidx[key][1],
        )

        # fig.update_layout(yaxis=)

        if title is not None:
            fig.update_layout(title=title)

    # fig['layout']['yaxis']['scaleanchor'] = 'x'
    # fig['layout']['yaxis']['scaleratio'] = 1

    fig.update_xaxes(dict(side='bottom', tickangle=-90, constrain='domain', scaleanchor='x', type='category'))
    fig.update_yaxes(dict(scaleanchor='x', scaleratio=1, type='category'))
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)')
    fig.show(renderer='firefox')
    return fig


def plotly_corrmat(df):
    """
    Plotly version of correlation matrix

    Parameters
    ----------
    df -

    Returns
    -------

    """
    xl = list(map(lambda x: f"{x[1]} {x[0]:>3}", df.columns.to_list()))
    yl = list(map(lambda x: f"{x[1]} {x[0]:>3}", df.index.to_list()))

    fig = px.imshow(df,
                    x=xl,
                    y=yl,
                    color_continuous_scale='rdbu_r',
                    color_continuous_midpoint=0.0)

    fig = go.Figure(data=go.Heatmap(
        z=df.values,
        x=xl,
        y=yl,
        colorscale='rdbu_r',
        zmin=-1.0,
        zmax=1.0,
    ))
    fig.update_layout(yaxis=dict(scaleanchor='x', scaleratio=1))
    fig.show(renderer='firefox')

    fig.update_coloraxes()

    fig.update_xaxes(
        tickangle=90,
        title_text="stimulus",
        title_font={"size": 12},
        title_standoff=25)

    fig.show(renderer='firefox')
    return fig


def plotly_corrscatter(df):
    """

    Parameters
    ----------
    df : long datafrme, with columns

    Returns
    -------

    """
    stimname = list(map(lambda x: f"{x.stimname_row:>10}, {x.stimname_col:>10}",
                        df.itertuples(index=False)))
    df['stimname'] = stimname

    stimfilt_list = [('1-6ol', '1-6ol'),
                     ('EP', 'EP'),
                     ('1-6ol+EP', '1-6ol+EP'),
                     ('1-6ol', 'EP'),
                     ('1-6ol', '1-6ol+EP'),
                     ('EP', '1-6ol+EP')]
    stimname_filt = list(map(lambda x: f"{x[0]:>10}, {x[1]:>10}", stimfilt_list))
    fig = px.violin(df.loc[df['stimname'].isin(stimname_filt)],
                    y="value",
                    x= 'conc', #"movie_conc",
                    box=True,
                    points="all",
                    color='conc',
                    facet_col='stimname',
                    violinmode='overlay')
    fig.update_traces(meanline_visible=True,
                      width=1,
                      bandwidth=.05,
                      points='all',
                      side='positive',
                      pointpos=-0.2,
                      jitter=0.05,
                      box=dict(line=dict(color='black')),
                      marker=dict(opacity=0.7, line=dict(width=0)))
    fig.update_xaxes(showgrid=True)
    fig.update_yaxes(showgrid=True)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

    # fig.update_layout(
    #    title=tstr)
    fig.show(renderer='firefox')
    return fig


# def total_fluorescence(df_frametimeinfo, totF):

def dfof_facet_plot(dfof, df_frametimeinfo, cids=None, yrange=(-0.5, 1)):
    if cids is None:
        cid0 = 20
        dcid = 40
        cids = list(range(cid0, cid0 + dcid))

    df = pd.DataFrame(dfof.T, index=df_frametimeinfo['frame_times'])
    fig = px.line(df,
                  x=df.index,
                  y=cids,
                  facet_row='variable',
                  facet_col=df_frametimeinfo['scope_pulse_idx'],
                  facet_row_spacing=0.03,  # default is 0.07 when facet_col_wrap is used
                  facet_col_spacing=0.03,  # default is 0.03,
                  labels=dict(x="time (s)", y=None),
                  width=600, height=60 * len(cids) + 50
                  )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.for_each_trace(lambda t: t.update(name=""))
    fig.update_yaxes(title_text=None, tick0=0.0, dtick=2, ticks='outside', col=1)
    fig.update_xaxes(matches=None)
    fig.update_yaxes(range=yrange)
    fig.update_layout(margin=dict(autoexpand=False, l=40, r=20, b=40, t=20, pad=0))
    # fig.show(renderer='firefox')
    return fig


def dfof_odor_facet(dfof, df_frametimeinfo, df_stimulus, cids=None, yrange=(-0.5, 2)):
    df = df_frametimeinfo.join(df_stimulus.set_index('stim_idx'), on='trialtrace_idx', how='inner')
    df['time0'] = df['frame_times'] - df['stim_ict']
    idx = df.index
    df.reset_index(inplace=True)

    if cids is None:
        cids = list(range(10))
    df_dfof = pd.DataFrame(dfof[:, idx].T)

    fig = px.line(df_dfof,
                  x=df['time0'],
                  y=cids,
                  facet_row='variable',
                  facet_col=df['stimname'],
                  color=df['trialtrace_idx'],
                  facet_row_spacing=0.01,  # default is 0.07 when facet_col_wrap is used
                  facet_col_spacing=0.01,  # default is 0.03,
                  labels=dict(x="time (s)", y=None),
                  width=400, height=100 * len(cids) + 50
                  )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.for_each_trace(lambda t: t.update(name=""))
    fig.update_yaxes(title_text=None, tick0=0.0, dtick=2, ticks='outside', col=1)
    fig.update_xaxes(matches=None)
    fig.update_yaxes(range=yrange)
    fig.update_layout(margin=dict(autoexpand=False, l=40, r=40, b=40, t=20, pad=0))
    # fig.show(renderer='firefox')
    return fig
