import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pathlib import Path
import pandas as pd

import pprint
from drosolf import orns, pns
import scipy.io as sio
import importlib
import textwrap
from rastermap.mapping import Rastermap
from itertools import product
from rearing import hallem, io, meta

pp = pprint.PrettyPrinter(indent=2, compact=True)

# make instance of HallemSet, for easy lookup of hi, odor names, and abbreviations
Hallem = hallem.HallemSet()

# %% Set dataset info and directory locations
# set dataset info

expt = dict()
expt['date'] = '2021-07-27'
expt['fly'] = 1
expt['microscope'] = 'galvo-res'

# make directory paths
PRJ_DIR = Path.cwd()
FLY_DIR = PRJ_DIR.joinpath('data', 'processed_data', expt['date'], str(expt['fly']))

# load metadata from yaml file in fly folder
fly_id = f"fly{expt['fly']}"
panel_idx = 0
# %%
metadata = io.load_yaml(FLY_DIR.parent.joinpath("metadata.yaml"))
acquisitions = metadata['Acquisitions'][expt['fly']]
odor_panel = metadata['Olfactometer']['pin odors']
#%%
expt['rearing'] = metadata['Fly info'][expt['fly']]['reared']

# session info
df_sessioninfo = meta.get_session_info(acquisitions, expt)
df_sessioninfo.rename_axis('iacq')
df_sessioninfo.insert(2, 'rearing', expt['rearing'])


if not FLY_DIR.is_dir():
    FLY_DIR.mkdir()

with pd.ExcelWriter(FLY_DIR.joinpath('flyinfo.xlsx')) as writer:
    df_sessioninfo.to_excel(writer, sheet_name='flyinfo', index=False)
#
# olfactometer info
df_olfactometer = meta.get_olfactometer_info(odor_panel)
hi_list = set(df_olfactometer.hi.to_list())
df_stimtype = meta.get_stimtype(hi_list)
#%%

# save metadata.xlsx
print(f"\n_________________DATASET_________________")
print(f"date: {expt['date']}\n"
      f"fly: {expt['fly']}\n"
      f"microscope: {expt['microscope']}")
print(f"\nDF_SESSIONINFO _____________")
pp.pprint(df_sessioninfo)

print(f"DF_OLFACTOMETER _____________")
pp.pprint(df_olfactometer)

# iterate through movies
for iacq in range(len(acquisitions)):
    expt['movie'] = acquisitions[iacq]['thorimage']  # name of subfolder in FLY_DIR
    print(expt['movie'])
    PROC_DIR = FLY_DIR.joinpath(expt['movie'])
    S2P_DIR = PROC_DIR.joinpath('suite2p', 'plane0')

    if not PROC_DIR.exists():
        PROC_DIR.mkdir()
    # print(f"\n{expt['movie']} : DF_STIMULUS")

    stimulus = acquisitions[iacq]['stimulus']
    if stimulus is not None:
        df_stimulus0 = pd.DataFrame({'pinA': stimulus['channelA'], 'pinB': stimulus['channelB']})
        n_stim = df_stimulus0.shape[0]
        # pp.pprint(df_stimulus0)

        df_chanA = meta.get_channel_info(df_stimulus0.pinA.to_list(), df_olfactometer)
        df_chanB = meta.get_channel_info(df_stimulus0.pinB.to_list(), df_olfactometer)
        df_stimulus = df_chanA.join(df_chanB, lsuffix='_a', rsuffix='_b')
        df_stimulus = df_stimulus.join(df_stimtype.set_index(['hi1', 'hi2']), on=['hi_a', 'hi_b'])
        df_stimulus.insert(0, 'stim_idx', np.arange(n_stim)+1)

        do_save = True

        if do_save:
            with pd.ExcelWriter(PROC_DIR.joinpath('metadata.xlsx')) as writer:
                df_sessioninfo.iloc[[iacq], :].to_excel(writer, sheet_name='sessioninfo')
                df_olfactometer.to_excel(writer, sheet_name='olfactometer', index=False)
                df_stimulus0.to_excel(writer, sheet_name='stimulus0', index=False)
                df_stimulus.to_excel(writer, sheet_name='stimulus', index=False)
                df_stimtype.to_excel(writer, sheet_name='stimtype', index=False)
