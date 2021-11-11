import argparse
import itertools
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats, signal
from scipy.spatial.distance import pdist, squareform
from sklearn import preprocessing
from sklearn import svm
from upsetplot import UpSet

from dataload import SparsityAnalysis, load_all_flydata, make_fly_list
from rearing import hallem, io, corrmat
from rearing.ca import respvec, trace, vis
import upsetplot

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'legend.fontsize': 8,
                     'axes.labelsize': 8,
                     'axes.titlesize': 10,
                     'xtick.labelsize': 8,
                     'ytick.labelsize': 8,
                     'font.family': ['sans-serif'],
                     'figure.dpi': 400.0
                     })

Hallem = hallem.HallemSet()
print("IT'S WORKING")


# ----------------------------
# Load db.csv
# ----------------------------

def main():
    save_plots = True

    db = io.load_db()
    fly_list = make_fly_list(db)

    for fly in fly_list:
        for mov in fly.expt_list:
            subdirlist = [item.name for item in mov.PROC_DIR().glob('__*')]

            for loadsubdir in subdirlist:
                # -----------------------------------------------
                # loop through parameter subfolders
                # -----------------------------------------------

                # load yaml options
                dfof_ops = np.load(mov.PROC_DIR().joinpath(loadsubdir, "dfof_ops.npy"),
                                   allow_pickle=True).item()

                # load processed data
                procdat = io.load_processed_data(mov.PROC_DIR().joinpath(loadsubdir))

                # title string
                tstr = f"date: {fly.date}, fly: {fly.fly_num}, movies: {mov.movie} [{procdat['params']['movie_conc']}] " \
                       f"\nrearing: {fly.metadata['Fly info'][fly.fly_num]['reared']}" \
                       f"\nwin={dfof_ops['win']}, prctile={dfof_ops['percentile']}"

                df_binresp = respvec.get_binresponse(procdat['responses']['df_peakresp'],
                                                     procdat['responses']['baseline_std'],
                                                     n_std=3, thr=0.5)
                df_tuning = df_binresp.groupby(by='stimname', axis=1).any()

                df_sets = df_tuning.drop(columns=['pfo', '1-5ol']).value_counts().sort_index()
                fig_upset = upsetplot.plot(df_sets, show_percentages=True,
                                           show_counts=True,
                                           orientation='horizontal')
                fig = fig_upset['matrix'].get_figure()
                fig.set_size_inches(12, 7)
                plt.suptitle(tstr)
                plt.show()

                if save_plots:
                    savenam = f"UpSetPlot__win{dfof_ops['win']}__{dfof_ops['percentile']}.pdf"
                    fig.savefig(mov.PROC_DIR().joinpath(loadsubdir, savenam))

# datadir = args.datadir
# n_flies = args.n_flies
# t_experiment = args.t_experiment  # all in seconds
# t_motor_on = args.t_motor_on
# motor_duration = args.motor_duration
# t_led_on = args.t_led_on
# led_duration = args.led_duration
# write_data = not args.selfcheck
