from pathlib import Path
from collections import namedtuple
from dataclasses import dataclass, field
from typing import List
from pathlib import Path
import pandas as pd
import numpy as np
import suite2p
from rearing import hallem, io
from rearing.ca import respvec, trace, vis
from scipy import stats, signal
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import svm
from sklearn import preprocessing
import math
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from sklearn.svm import SVC
from scipy.signal import convolve
from itertools import compress
import copy

Hallem = hallem.HallemSet()


@dataclass
class SparsityAnalysis:
    sparsity: np.ndarray = field(metadata={'axis0': 'len(n_std)', 'axis1': '# of trials/stimuli'})
    n_std: np.ndarray
    stim_idx: List
    trial_dict: dict = field(default_factory=dict)
    df: pd.DataFrame = field(init=False)
    df_melted: pd.DataFrame = field(init=False)

    def __post_init__(self):
        self.df = pd.DataFrame(self.sparsity, index=self.n_std, columns=self.stim_idx)
        self.df.index.name = 'n_std'
        self.trial_dict.update({'stim_idx': self.stim_idx})

        if len(self.trial_dict) > 0:
            mux = pd.DataFrame(self.trial_dict)
            self.df.columns = pd.MultiIndex.from_frame(mux)

            self.df_melted = pd.melt(self.df.reset_index(), id_vars='n_std', value_name='sparsity')


@dataclass
class Fly:
    date: str
    fly_num: int
    rearing: str
    movies: List[str] = field(default_factory=list, compare=False)
    expt_list: List = field(default_factory=list, compare=False)
    PRJ_DIR: Path = field(default_factory=Path.cwd, compare=False)
    FLY_DIR: Path = field(init=False, compare=False)
    metadata: dict = field(default_factory=dict, compare=False)

    def __post_init__(self):
        self.FLY_DIR = self.PRJ_DIR.joinpath('data', 'processed_data', self.date, str(self.fly_num))
        self.load_metadata_yaml()
        self.expt_list = [None] * len(self.movies)
        for idx, mov in enumerate(self.movies):
            expt = Experiment(self, mov)
            self.expt_list[idx] = expt

    def load_metadata_yaml(self):
        """
        Load metadata.yaml from fly folder


        :return: metadata, contains all data from "metadata.yaml" file.
        :rtype: dict
        """
        try:
            metadata = io.load_yaml(self.FLY_DIR.joinpath("metadata.yaml"))
            self.metadata = metadata
        except (ValueError, KeyError, OSError,
                RuntimeError, TypeError, NameError):
            print(f"ERROR: Can't load metadata.yaml from {self.__repr__}")

    def __str__(self):
        str1 = f"Fly(date={self.date}, fly_num={self.fly_num}, rearing={self.rearing}, movies={self.movies})"
        str2 = f"\n\tPRJ_DIR = {self.PRJ_DIR}" + f"\n\tFLY_DIR = {self.FLY_DIR}"
        str3 = "\n".join([f"\t\t{k}:{{...}}" for k in self.metadata.keys()])
        return "\n" + str1 + str2 + "\n\tmetadata: { \n" + str3 + "\t}"


@dataclass
class Experiment:
    fly_parent: Fly
    movie: str

    # PROC_DIR: Path = field(init=False, compare=False)
    # S2P_DIR: Path = field(init=False, compare=False)

    def PROC_DIR(self) -> Path:
        return self.fly_parent.FLY_DIR.joinpath(self.movie)

    def S2P_DIR(self, combined=False):
        if combined:
            s2p_dir = self.PROC_DIR().joinpath('suite2p', 'combined')
            if s2p_dir.exists():
                return s2p_dir
        else:
            return self.PROC_DIR().joinpath('suite2p', 'plane0')
    def __hash__(self):
        return hash((self.date, self.fly_num, self.movie))

    def basefilename(self):
        return f"{self.fly_parent.date}__{self.fly_parent.fly_num:02}__{self.movie}"

    def __repr__(self):
        return f"Experiment(  movie={self.movie}, " \
               f"  fly_parent=Fly( date={self.fly_parent.date}, fly_num={self.fly_parent.fly_num})  )"

    def load_ti(self):
        # load timing/stimulus data
        ti = io.load_ti(self.PROC_DIR().joinpath('ti.mat'))
        return ti

    def load_df_stimulus(self):
        df_stimulus = pd.read_excel(self.PROC_DIR().joinpath('metadata.xlsx'), sheet_name='stimulus')
        conc = get_unique_except(df_stimulus['conc_a'], df_stimulus['conc_b'])
        df_stimulus['conc'] = conc
        return df_stimulus

    def load_suite2p_data(self):
        dat0 = io.load_s2p_files(self.S2P_DIR().joinpath('stat.npy'))
        dat0 = io.preprocess_s2p_dat(dat0)
        return dat0

    def load_frametimeinfo(self):
        df_frametimeinfo = pd.read_csv(self.PROC_DIR().joinpath('frametimeinfo.csv'))
        return df_frametimeinfo


@dataclass(init=True, repr=True, eq=True, order=True, unsafe_hash=False, frozen=False)
class MultiChanOdorStim:
    hi1: List[int] = field(default_factory=list, compare=True)
    hi2: List[int] = field(default_factory=list, compare=True)


@dataclass(init=True, repr=True, eq=True, order=True, unsafe_hash=False, frozen=False)
class OdorStim:
    hi: int
    conc: int
    odor: str = field(init=False, compare=False)
    abbrev: str = field(init=False, compare=False)

    def __post_init__(self):
        self.odor = Hallem.hi2odor(self.hi)
        self.abbrev = Hallem.hi2abbrev(self.hi)

    def __str__(self):
        return "{:<12}[{:>3}]".format(self.abbrev, self.conc)


# load all fly data into desired structure

def get_unique_except(a, b, except_val=0):
    x = [list(set(item)) for item in zip(a, b)]
    for i, val in enumerate(x):
        if len(val) > 1:
            temp = [item for item in val if item != except_val]
            x[i] = temp
    for \
            i, val in enumerate(x):
        if len(val) == 1:
            x[i] = val[0]
    return x


def make_fly_list(db0):
    flies = []
    for row in db0.itertuples():
        print(row)
        movies = list(map(str.strip, row.movies.split(",")))
        fly = Fly(row.date, row.fly_num, row.rearing, movies)
        print(fly)

        expt_list = []
        for mov in movies:
            # Experiment instance for each movie
            expt = Experiment(fly, mov)
            expt_list.append(expt)

        setattr(fly, 'expt_list', expt_list)
        flies.append(fly)
    return flies


def load_all_flydata(db0, load_s2p_dat=True):
    for row in db0.itertuples():
        print(row)
        movies = list(map(str.strip, row.movies.split(",")))
        fly = Fly(row.date, row.fly_num, row.rearing, movies)
        print(fly)

        expt_list = []
        for mov in movies:
            # Experiment instance for each movie
            expt = Experiment(fly, mov)

            # load timing/stimulus data
            ti = io.load_ti(expt.PROC_DIR().joinpath('ti.mat'))

            # load stimulus trial dataframe
            df_stimulus = pd.read_excel(expt.PROC_DIR().joinpath('metadata.xlsx'), sheet_name='stimulus')
            conc = get_unique_except(df_stimulus['conc_a'], df_stimulus['conc_b'])
            df_stimulus['conc'] = conc

            # load frame info
            df_frametimeinfo = pd.read_csv(expt.PROC_DIR().joinpath('frametimeinfo.csv'))
            if 'isbadframe' not in df_frametimeinfo.columns:
                df_frametimeinfo['isbadframe'] = 0

                temp = df_frametimeinfo.loc[
                    df_frametimeinfo['olf_pulse_idx'] > 0, ['olf_pulse_idx',
                                                            'scope_pulse_idx']].drop_duplicates().set_index(
                    'olf_pulse_idx')
                df_stimulus = df_stimulus.join(temp, on='stim_idx')
                df_stimulus.rename(columns={"scope_pulse_idx": "block_idx"}, inplace=True)

            fti = df_frametimeinfo.to_records()

            # load suite2p data from suite2p folder, attach to expt object
            # dat0 = Suite2pData(*suite2p.gui.io.load_files(expt.S2P_DIR().joinpath('stat.npy')))
            dat0 = io.load_s2p_files(expt.S2P_DIR().joinpath('stat.npy'))
            dat0 = io.preprocess_s2p_dat(dat0)

            setattr(expt, 'ti', ti)
            setattr(expt, 'df_stimulus', df_stimulus)
            setattr(expt, 'df_frametimeinfo', df_frametimeinfo)
            setattr(expt, 'fti', fti)

            if load_s2p_dat:
                setattr(expt, 's2p_dat', dat0)

            expt_list.append(expt)

        setattr(fly, 'expt_list', expt_list)
        flies.append(fly)
    return flies


def good_trials(df_frametimeinfo_):
    """
    Takes df_frametimeinfo (loaded from frametimeinfo.csv, should have column 'isbadframe'), and returns a
    DataFrame indicating which stimuli/trials (in df_stimulus) are good.

    :param df_frametimeinfo_: Loaded from frametimeinfo.csv, not cropped to exclude bad frames
    :type df_frametimeinfo_: pd.DataFrame
    :return: df_good_trials: columns = {"olf_pulse_idx", "scope_pulse_idx", "isgoodtrial"}
    :rtype:
    """
    df_frametimeinfo = df_frametimeinfo_.loc[df_frametimeinfo_['isbadframe'] == 0]

    # dataframe w/ columns "olf_pulse_idx", "scope_pulse_idx" for all trials in df_frametimeinfo_
    df_all_trials = df_frametimeinfo_.loc[df_frametimeinfo_['olf_pulse_idx'] > 0,
                                          ['olf_pulse_idx', 'scope_pulse_idx']].drop_duplicates()

    # dataframe w/ columns "olf_pulse_idx", "scope_pulse_idx", but w/ rows corresponding to "bad stimuli"
    # (occurred during "isbadframe" time points) dropped
    df_good_trials = df_frametimeinfo.loc[df_frametimeinfo['olf_pulse_idx'] > 0,
                                          ['olf_pulse_idx', 'scope_pulse_idx']].drop_duplicates()

    df_all_trials['isgoodtrial'] = df_all_trials.isin(df_good_trials).any(axis=1).to_list()
    df_all_trials.reset_index(drop=True, inplace=True)
    return df_all_trials


def s2p_frametrim(idx, s2p_dat_, return_changed_only=False):
    """
    Given a dict of suite2p data and frametime indices to keep, trim timeseries data
    (Fneu, Fcell, Spks, F, Fc) to drop bad frametimes.

    Trimmed timeseries can be returned by themselves, or in a dict with the other suite2p data (like input)

    :param idx: indices for good frames, correspond to df_frametimeinfo['isbadframe']==0
    :type idx: np.ndarray
    :param s2p_dat_: contains suite2p outputs
    :type s2p_dat_: dict
    :param return_changed_only:
    :type return_changed_only: bool

    :return: s2p_dat
    :rtype: dict

    """
    varnames = ['Fcell', 'Fneu', 'Spks', 'F', 'Fc']
    s2p_dat = {}
    for item in varnames:
        s2p_dat[item] = s2p_dat_[item][:, idx]
    if return_changed_only:
        return s2p_dat
    else:
        s2p_dat_ = s2p_dat_.copy()
        s2p_dat_.update(s2p_dat)
        return s2p_dat_


def adjust_for_bad_frames(df_frametimeinfo_, df_stimulus_, s2p_dat0=None):
    """
    Trims df_frametimeinfo and fluorescence traces in s2p_dat to exclude bad frames
    (based on "isbadframe" column). The index of df_frametimeinfo is reset to 0.

    Also adds column "isgoodtrial" to df_stimulus, and trims timeseries data in s2p_dat if passed.

    :param df_frametimeinfo_: original df_frametimeinfo, loaded from 'df_frametimeinfo.csv' (has column 'isbadframe')
    :type df_frametimeinfo_: pd.DataFrame
    :param df_stimulus_: loaded from metadata.xlsx, sheetname='stimulus'
    :type df_stimulus_: pd.DataFrame
    :param s2p_dat0: contains suite2p outputs (from .npy files),
    :type s2p_dat0: dict
    :return:
    :rtype:
    """
    if 'isbadframe' in df_frametimeinfo_.columns:
        # idx : which frames to keep
        idx = (df_frametimeinfo_['isbadframe'] == 0).to_numpy()

        # trim df_frametimeinfo to include only good frames
        df_frametimeinfo = df_frametimeinfo_.loc[df_frametimeinfo_['isbadframe'] == 0]
        df_frametimeinfo.reset_index(drop=True, inplace=True)

        # add "isgoodtrial" column to df_stimulus
        df_good_trials = good_trials(df_frametimeinfo_)
        df_good_trials.set_index(['olf_pulse_idx', 'scope_pulse_idx'], inplace=True)
        df_stimulus = df_stimulus_.join(df_good_trials, on=['stim_idx', 'block_idx'])

        # trim timeseries in s2p_dat if passed, then return
        if s2p_dat0 is None:
            return df_frametimeinfo, df_stimulus
        else:
            s2p_dat = s2p_frametrim(idx, s2p_dat0, return_changed_only=False)
            return df_frametimeinfo, df_stimulus, s2p_dat

    # no bad frames
    else:
        if s2p_dat0 is None:
            return df_frametimeinfo_, df_stimulus_
        else:
            return df_frametimeinfo_, df_stimulus_, s2p_dat0


def stim2cat(df_stimulus0):
    df_stimulus = copy.deepcopy(df_stimulus0)

    ord_all = ['empty', 'pfo', '1-5ol', '1-6ol', '1-5ol+1-6ol', 'EP']
    ord_mov = list(filter(lambda x: x in df_stimulus.stimname.unique(), ord_all))
    cat_type = pd.api.types.CategoricalDtype(categories=ord_mov, ordered=True)
    df_stimulus['stimname'] = df_stimulus['stimname'].astype(cat_type)
    return df_stimulus

# %%

# %%
# stim1 = OdorStim(73, -3)
# stim2 = OdorStim(73, -5)
# stim3 = OdorStim(73, -5)
# stim3.strtype = 'odor'
#
# proc_dir_rel = ('processed_data',)
#
# stimlist = [stim1, stim2, stim3]
