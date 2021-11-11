#  STEP 1. Setup

import copy
import json
import pprint as pp

from pathlib import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import plotly.express as px
import plotly.io as pio
from scipy import signal
from scipy.stats import zscore

from rearing import tidy
from ryim.dirutils import *
from rearing.ca import respvec, trace, vis, represp
import rasterplot

import pyarrow.feather as feather
import seaborn as sns

pio.renderers.default = "browser"


# %%


def pooled_corrscatter(files):
    DATA_DIR = Path(r'./data/processed_data')

    tidy_files_all = sorted(list(DATA_DIR.rglob(r"**/tidy_corr_all.feather")))
    tidy_files_resp = sorted(list(DATA_DIR.rglob(r"**/tidy_corr_resp.feather")))

    tidy_files = tidy_files_resp
    tidy_list = [None] * len(tidy_files)
    for i, f in enumerate(tidy_files_all):
        tidy_df = feather.read_feather(f)

        # get experiment info
        mov_dir = parse_parent_movie(f)
        flydate = parse_parent_date(f).name
        flynum = parse_parent_fly(f).name
        mov = parse_parent_movie(f).name

        # load session info, determine rearing type, load stimulus timing info from sync_meta.json
        session_info = pd.read_excel(mov_dir.joinpath('metadata.xlsx'), sheet_name='sessioninfo').squeeze()
        rear_type = session_info['rearing']

        # set rearing condition column
        if 'pfo' in rear_type:
            tidy_df['rearing'] = 'pfo'
        else:
            tidy_df['rearing'] = 'odor'

        tidy_list[i] = tidy_df

    tidy_cat = pd.concat(tidy_list, axis=0)

    order = tidy_cat['stimname'].unique().tolist()
    # tidy_cat['stimname'] = list(zip(tidy_cat['stimname_row'], tidy_cat['stimname_col']))
    tidy_cat['stimname'] = list(map(lambda x: f"{x.stimname_row}, {x.stimname_col}",
                                    tidy_cat.itertuples(index=False)))

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)

    g = sns.catplot(data=tidy_cat.loc[tidy_cat['corrtype'] == 'cosine', :],
                    x="stimname", y="value", hue='rearing',
                    row='conc')

    sns.stripplot(data=tidy_cat.loc[tidy_cat['corrtype'] == 'cosine', :],
                  x="stimname", y="value", hue='rearing',
                  ax=ax,
                  dodge=True, order=order, palette="pastel")

    ax.set_ylim((-0.2, 1.1))
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8, rotation=45, ha="right",
                       rotation_mode="anchor")
    sns.despine()
    plt.show()
    return fig, ax
