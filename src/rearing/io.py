# import seaborn as sns
import re
import pandas as pd
import yaml
import suite2p.gui
import numpy as np
import scipy.io as sio
from pathlib import Path


def load_ti(fname_):
    ti0 = sio.loadmat(fname_)
    ti = {key: value.squeeze().tolist() for key, value in ti0.items() if isinstance(value, np.ndarray)}
    return ti


def load_yaml(fname):
    """
    Loads yaml from file.

    :param fname: yaml filepath
    :type fname: Path
    :return: meta - dictionary containing all the information in the yaml file
    :rtype: dict
    # """
    with open(fname, 'r') as yamlfile:
        meta = list(yaml.safe_load_all(yamlfile))
    meta = meta[0]
    # with open(fname, 'r') as stream:
    #     meta = yaml.load(stream, Loader=yaml.SafeLoader)
    return meta


def load_s2p_files(fname):
    """
    Loads suite2p data into dictionary, doesn't do anything w/ icell variable
    :param fname: filepath for stat.npy file
    :type fname: pathlib.PurePath
    :return: s2p_dat - dict w/ keys
                            ('stat', 'ops', 'Fcell', 'Fneu','Spks',
                              'iscell', 'probcell', 'redcell', 'probredcell', 'hasred')
    """
    vals = suite2p.gui.io.load_files(fname)
    keys = ('stat', 'ops', 'Fcell', 'Fneu', 'Spks', 'iscell', 'probcell', 'redcell', 'probredcell', 'hasred')
    s2p_dat = dict(zip(keys, vals))
    return s2p_dat


def load_db():
    db = pd.read_csv(Path.cwd().joinpath('src', 'db.csv'))
    return db


def load_processed_data(folder, as_dict=True):
    params = np.load(folder.joinpath('params.npy'), allow_pickle=True).item()
    dfof = np.load(folder.joinpath('dfof.npy'))

    dist = np.load(folder.joinpath('distances.npy'), allow_pickle=True).item()
    dist_triu = np.load(folder.joinpath('dist_triu.npy'), allow_pickle=True).item()
    dist_triu_asc = np.load(folder.joinpath('dist_triu.npy'), allow_pickle=True).item()
    responses = np.load(folder.joinpath('responses.npy'), allow_pickle=True).item()
    sparsity = np.load(folder.joinpath('sparsity.npy'), allow_pickle=True).item()

    dfof = np.load(folder.joinpath('dfof.npy'))
    dfof = np.load(folder.joinpath('dfof.npy'))

    if as_dict:
        return dict(params=params, dfof=dfof, sparsity=sparsity, responses=responses, dist=dist)
    else:
        #return params, dfof, sparsity, responses, dist, dist_triu, dist_triu_asc
        return params, dfof, sparsity, responses, dist


def preprocess_s2p_dat(data_dict, neu_coeff=0.7):
    iscell = data_dict['iscell']
    F = data_dict['Fcell'][iscell, :]
    Fneu = data_dict['Fneu'][iscell, :]
    Fc = F - neu_coeff * Fneu
    dat = data_dict.copy()
    dat.update({'iscell': iscell, 'F': F, 'Fneu': Fneu, 'Fc': Fc})
    return dat

# def parse_pin_lookup(pin_lookup_dict):
#     """
#     Takes a dictionary (usually from metadata.yaml file) of pin lookup info, and returns dataframe w/ odor
#     and concentration columns.
#     :param dict pin_lookup_dict: dictionary containing stimulus pins and text lookup info
#     """
#     df_pins = pd.DataFrame.from_dict(pin_lookup_dict, orient='index', columns=['txt'])
#     df_pins['odor'] = ''
#     df_pins['conc_log_10'] = ''
#     for index, value in df_pins['txt'].items():
#         if value is not None:
#             odor = value.split('(')[0].strip()
#             try:
#                 conc_log_10 = float(re.search(r"\((.*?)\)", value).group(1))
#             except:
#                 conc_log_10 = None
#             df_pins.at[index, 'odor'] = odor
#             df_pins.at[index, 'conc_log_10'] = conc_log_10
#         else:
#             df_pins.at[index, 'odor'] = None
#             df_pins.at[index, 'conc_log_10'] = None
#     df_pins.drop(labels='txt', axis=1, inplace=True)
#     return df_pins
#
#
# def pinlist_to_stimdf(df_pins, pins):
#     """
#     Given a dataframe with pin lookup info (index = arduino pin, columns = ['odor', 'conc_log_10']) and a list of pins
#     used in the arduino script for stimulus delivery, returns a dataframe  with columns = ['odor', 'conc_log_10'] for
#     every stimulus row.
#
#     :param df_pins: Dataframe with index = arduino pin, and columns ['odor', 'conc_log_10']; generated from io.parse_pin_lookup
#     :type df_pins: Pandas.Dataframe
#     :param pins: list of pins in order of stimulus delivery
#     :type pins: list
#     """
#
#     df = df_pins.loc[pins]
#     df['pins'] = pins
#     df = df[['pins', 'odor', 'conc_log_10']]
#     return df
