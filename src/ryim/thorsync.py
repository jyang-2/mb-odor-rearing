import numpy as np
import collections
import glob
from pathlib import Path
import json


clockrate = 20000000


#
# with h5py.File(file, 'r') as f:
#     scopePin = f['AI']['scopePin'][].flatten()
#     scopePin = (scopePin > 4)

def load_sync_meta(json_filepath: Path):
    """
    Load sync_meta.json into dict w/ fields:

        - block_fct : start time for each acquisition
        -  block_ict: end time for acquisition
        -   n_blocks: # of acquisitions
        -     n_stim: total # of stimuli presented
        -   stim_ict: stimulus onset times
        -   stim_fct: stimulus offset times

    Note: this is always loaded from the movie folder, not the block folder.

    If you are analyzing a block, then go to the parent folder to load, and note
    that it will include stimuli info for different blocks as well.

    """

    # load into dict ti (trial info)
    with open(json_filepath) as file:
        ti = json.load(file)

    # make everything into a 1d numpy array
    for k, v in ti.items():
        if isinstance(v, list):
            ti[k] = np.array(v)
    return ti


def get_xml_files(folder):
    thorsync_xml_files = sorted(list(folder.glob('**/Episode*.h5')))
    return thorsync_xml_files


def read_sync_data():
    sync_data = []
    sync_dt = []
    for epi_files in self.sync_episodes:
        for episode in epi_files:
            sync_data.append({})
            sync_dt.append({})
            print(episode)
            h5 = tables.open_file(episode)
            for el in h5.root.DI:
                sync_data[-1][el.name] = np.squeeze(el)
                sync_dt[-1][el.name] = self._find_dt(
                    el.name, len(sync_dt) - 1)
            h5.close()

    return sync_data, sync_dt
