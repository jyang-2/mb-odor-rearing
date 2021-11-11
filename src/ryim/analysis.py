from __future__ import annotations
from pathlib import Path
import pathlib
import pandas as pd
import numpy as np
import re
from collections import namedtuple
from abc import ABC, abstractmethod
from ryim import dirutils, thorsync
import pprint as pp
from natsort import natsorted, realsorted
import yaml
import json
from rearing.ca import respvec, trace, vis, represp

ocat = pd.CategoricalDtype(categories=['pfo', '1-5ol', 'VA', '1-6ol', 'EP', '1-6ol+EP'], ordered=True)


# useful functions
def all_plane_folders(suite2p_dir: Path | str) -> list(Path):
	"""
	Returns list of all plane folders in suite2p (name is of the form "plane*")

	:param suite2p_dir: parent suite2p folder
	"""
	suite2p_dir = Path(suite2p_dir)

	plane_folders = [item for item in suite2p_dir.glob('plane*') if item.is_dir()]
	plane_folders = realsorted(plane_folders, key=lambda x: x.name)
	return plane_folders


def plane_number(folder: Path | str) -> int | list[int]:
	"""
	Returns plane # from plane folder filepath.

	Assumes if folder is a string, it's just the folder name.

	:param folder:
	:return:
		int - plane number
	"""

	if isinstance(folder, Path):
		nums = re.findall(r'\d+', folder.name)
	elif isinstance(folder, str):
		nums = re.findall(r'\d+', folder)
	elif isinstance(folder, list):
		return [plane_number(item) for item in folder]
	return int(nums[0])


def postprocess_df_frametimeinfo(_df_frametimeinfo):
	"""Drop rows if multiplane"""
	if 'plane' in _df_frametimeinfo.columns:
		_df_frametimeinfo = _df_frametimeinfo.loc[_df_frametimeinfo['plane'] == 0]
		_df_frametimeinfo.reset_index(drop=True, inplace=True)  #
	return _df_frametimeinfo


def get_extra_metadata_info(_nt_meta_data, ocat=ocat):
	"""
	Returns frame rate, block concentration, and rearing condition

	:param _nt_meta_data: namedtuple containing dataset analysis metadata
	:param ocat: pd.CategoricalDtype, for ordering df_stimulus['stimname']
	:return: tuple(float, int|float, str)
	          fs - float, frame rate, computed from nt_meta_data['timestamps']
         blk_conc- int | float, concentration for movie, computed from df_stimulus['stimname']
		rear_type-  str, rearing condition (listed in metadata.xlsx, sheet=sessioninfo)

	"""

	_ti = _nt_meta_data.sync_meta
	_df_stimulus = _nt_meta_data.df_stimulus
	_df_frametimeinfo = _nt_meta_data.df_frametimeinfo

	_df_stimulus['blk_conc'] = min(_df_stimulus.loc[:, ['conc_a', 'conc_b']].mode(axis=0).min(axis=0).tolist())
	_df_stimulus['stim_ict'] = _ti['stim_ict'][_df_stimulus['stim_idx'].to_numpy() - 1]
	_df_stimulus['scope_pulse_idx'] = np.interp(_df_stimulus['stim_ict'],
	                                            _df_frametimeinfo['frame_times'],
	                                            _df_frametimeinfo['scope_pulse_idx']).astype('int')
	_df_stimulus['stimname'] = _df_stimulus['stimname'].astype(ocat)

	_df_stimulus['reps'] = ''
	for ustim in _df_stimulus['stimname'].unique():
		mask = _df_stimulus['stimname'] == ustim
		_df_stimulus.loc[mask, 'reps'] = np.arange(np.count_nonzero(mask * 1))

	_session_info = {key: value for key, value in _nt_meta_data.session_info.to_dict().items() if key.isidentifier()}

	_fs = trace.compute_fs(_nt_meta_data.timestamps)
	_blk_conc = _df_stimulus.blk_conc.unique()[0]
	_rear_type = _session_info['rearing']

	return _fs, _blk_conc, _rear_type


class MetaDataProcessor(ABC):
	"""
	Abstract class for fixing and adding data to metadata files after loading w/ an instance of MetaDataLoader.
    """

	ProcessedMetaData = namedtuple("MetaData",
	                               ["df_stimulus",
	                                "df_frametimeinfo",
	                                "timestamps",
	                                "session_info",
	                                "xml_meta",
	                                "sync_meta",
	                                "meta_yaml",
	                                "fs",
	                                "blk_conc",
	                                "rear_type"
	                                ])

	@classmethod
	def process_metadata(cls, _nt_meta_data):
		"""
        Class method for adding information to nt_meta_data returned by an instance of MetaDataLoader
        :param _nt_meta_data:
        :return: ProcessedMetaData: named tuple w/ added fields
                        fs : frame rate
                        blk_conc : concentration of block
                        rear_type : rearing conditions
        """
		_ti = _nt_meta_data.sync_meta
		_df_stimulus = _nt_meta_data.df_stimulus
		_df_frametimeinfo = _nt_meta_data.df_frametimeinfo

		_df_stimulus['blk_conc'] = min(_df_stimulus.loc[:, ['conc_a', 'conc_b']].mode(axis=0).min(axis=0).tolist())
		_df_stimulus['stim_ict'] = _ti['stim_ict'][_df_stimulus['stim_idx'].to_numpy() - 1]
		_df_stimulus['scope_pulse_idx'] = np.interp(_df_stimulus['stim_ict'],
		                                            _df_frametimeinfo['frame_times'],
		                                            _df_frametimeinfo['scope_pulse_idx']).astype('int')
		_df_stimulus['stimname'] = _df_stimulus['stimname'].astype(ocat)

		_df_stimulus['reps'] = ''
		for ustim in _df_stimulus['stimname'].unique():
			mask = _df_stimulus['stimname'] == ustim
			_df_stimulus.loc[mask, 'reps'] = np.arange(np.count_nonzero(mask * 1))

		_fs = trace.compute_fs(_nt_meta_data.timestamps)
		_blk_conc = _df_stimulus.blk_conc.unique()[0]

		_session_info = _nt_meta_data['session_info']
		_rear_type = _session_info['rearing']

		procmeta: MetaDataProcessor.MetaData = cls.ProcessedMetaData(_df_stimulus,
		                                                             _df_frametimeinfo,
		                                                             _nt_meta_data.timestamps,
		                                                             _session_info,
		                                                             _nt_meta_data.xml_meta,
		                                                             _nt_meta_data.sync_meta,
		                                                             _nt_meta_data.meta_yaml,
		                                                             _fs,
		                                                             _blk_conc,
		                                                             _rear_type)
		return procmeta


class MetaDataLoader(ABC):
	"""
	Abstract class for loading metadata files (implementation depends on multiplane vs single, block vs movie
	"""
	MetaData = namedtuple("MetaData", ["df_stimulus", "df_frametimeinfo", "timestamps",
	                                   "session_info", "xml_meta", "sync_meta", "meta_yaml"])

	def __init__(self, suite2p_dir, mov_dir, fly_dir, date_dir):
		self.suite2p_dir = suite2p_dir
		self.mov_dir = mov_dir
		self.fly_dir = fly_dir
		self.date_dir = date_dir

	@classmethod
	@abstractmethod
	def init_from_statfile(cls, statfile):
		"""

		:param statfile: filepath to something in the same folder as suite2p/**/stat.npy
		"""
		suite2p_dir = dirutils.parse_parent_suite2p(statfile)
		blk_dir = dirutils.parse_parent_block(statfile)
		mov_dir = dirutils.parse_parent_movie(statfile)
		fly_dir = dirutils.parse_parent_fly(statfile)
		date_dir = dirutils.parse_parent_date(statfile)
		return cls(suite2p_dir, mov_dir, fly_dir, date_dir)

	def load_metadata_yaml(self):
		"""loads metadata.yaml file from date directory"""
		yaml_file = self.date_dir.joinpath("metadata.yaml")

		with open(yaml_file, 'r') as f:
			metadata_yaml = list(yaml.safe_load_all(f))
		return metadata_yaml[0]

	def load_session_info(self):
		""" Rearing info should be here """
		metadata_xlsx_filepath = self.mov_dir.joinpath('metadata.xlsx')
		session_info = pd.read_excel(metadata_xlsx_filepath, sheet_name='sessioninfo').squeeze()
		return session_info

	def load_sync_meta(self):
		"""Load thorsync-based info from sync_meta.json file"""
		return thorsync.load_sync_meta(self.mov_dir.joinpath("sync_meta.json"))

	def load_xml_meta(self):
		"""Load movie info from xml_meta.json"""
		xml_meta_filepath = self.mov_dir.joinpath("xml_meta.json")
		with open(xml_meta_filepath) as file:
			xml_meta_ = json.load(file)
		return xml_meta_

	@abstractmethod
	def load_df_frametimeinfo(self):
		"""Load df_frametimeinfo.csv"""
		pass

	@abstractmethod
	def load_df_stimulus(self):
		"""Load df_stimulus from metadata.xls or df_stimulus.csv"""
		pass

	@abstractmethod
	def load_timestamps(self):
		"""Load timestamps from timestamps.npy"""
		pass

	def load(self):
		"""Load all metadata"""
		_meta_yaml = self.load_metadata_yaml()
		_sync_meta = self.load_sync_meta()
		_xml_meta = self.load_xml_meta()
		_session_info = self.load_session_info()
		_df_stimulus = self.load_df_stimulus()
		_df_frametimeinfo = self.load_df_frametimeinfo()
		_timestamps = self.load_timestamps()

		if 'plane' in _df_frametimeinfo.columns:
			_df_frametimeinfo = _df_frametimeinfo.loc[_df_frametimeinfo['plane'] == 0]
			_df_frametimeinfo.reset_index(drop=True, inplace=True)  #
		return self.MetaData(_df_stimulus, _df_frametimeinfo, _timestamps, _session_info, _xml_meta, _sync_meta,
		                     _meta_yaml)


class BlockMetaDataLoader(MetaDataLoader):
	"""
	Data loader class for loading from movie blocks.

	"""

	def __init__(self, suite2p_dir, blk_dir, mov_dir, fly_dir, date_dir):
		super().__init__(suite2p_dir, mov_dir, fly_dir, date_dir)
		self.blk_dir = blk_dir

	@classmethod
	def init_from_statfile(cls, statfile):
		suite2p_dir = dirutils.parse_parent_suite2p(statfile)
		blk_dir = dirutils.parse_parent_block(statfile)
		mov_dir = dirutils.parse_parent_movie(statfile)
		fly_dir = dirutils.parse_parent_fly(statfile)
		date_dir = dirutils.parse_parent_date(statfile)
		return cls(suite2p_dir, blk_dir, mov_dir, fly_dir, date_dir)

	def load_df_frametimeinfo(self):
		"""Loads <mov_dir>/<blk_dir>/frametimeinfo.csv"""
		_df_frametimeinfo = pd.read_csv(self.blk_dir.joinpath('frametimeinfo.csv'))
		return postprocess_df_frametimeinfo(_df_frametimeinfo)

	def load_df_stimulus(self):
		return pd.read_csv(self.blk_dir.joinpath('metadata.csv'))

	def load_timestamps(self):
		return np.load(self.blk_dir.joinpath('timestamps.npy'))

	def load(self):
		return super().load()


class MovieMetaDataLoader(MetaDataLoader):
	"""
	Data loader class for loading from movies.

	"""
	def __init__(self, suite2p_dir, mov_dir, fly_dir, date_dir):
		super().__init__(suite2p_dir, mov_dir, fly_dir, date_dir)

	@classmethod
	def init_from_statfile(cls, statfile):
		return super().init_from_statfile(statfile)

	def load_df_frametimeinfo(self):
		"""Loads <mov_dir>/frametimeinfo.csv"""
		_df_frametimeinfo = pd.read_csv(self.mov_dir.joinpath('frametimeinfo.csv'))
		return postprocess_df_frametimeinfo(_df_frametimeinfo)

	def load_df_stimulus(self):
		return pd.read_excel(self.mov_dir.joinpath('metadata.xlsx'), sheet_name='stimulus')

	def load_timestamps(self):
		return np.load(self.mov_dir.joinpath('timestamps.npy'))

	def load(self):
		return super().load()


class MovieMetaDataLoader_1(MovieMetaDataLoader):
	MetaData = namedtuple("MetaData", ["df_stimulus", "df_frametimeinfo", "timestamps"])

	def load(self):
		_df_stimulus = self.load_df_stimulus()
		_df_frametimeinfo = self.load_df_frametimeinfo()
		_timestamps = self.load_timestamps()

		if 'plane' in _df_frametimeinfo.columns:
			_df_frametimeinfo = _df_frametimeinfo.loc[_df_frametimeinfo['plane'] == 0]
			_df_frametimeinfo.reset_index(drop=True, inplace=True)  #
		return self.MetaData(_df_stimulus, _df_frametimeinfo, _timestamps)


# %% Suite2p Data loaders
class Suite2pLoader:
	"""
	Parent class, basic Suite2p data loader
	"""
	Suite2pOutputs = namedtuple("Suite2pOutputs", ["F", "Fneu", "stat", "ops", "iscell"])

	@staticmethod
	def statfile_parent(statfile: str | Path):
		"""
		:param statfile: Union[str, Path] statfile: filepath to a "**/stat.npy" file or similar
		:return:
			- Path: parent folder from which to load suite2p results
		"""
		if isinstance(statfile, str):
			statfile = Path(statfile)
		return statfile.parent

	@staticmethod
	def load_from_folder(folder):
		"""
		Loads suite2p output files from a single folder.

		"""
		stat = np.load(folder.joinpath('stat.npy'), allow_pickle=True)
		ops = np.load(folder.joinpath('ops.npy'), allow_pickle=True).item()
		F = np.load(folder.joinpath('F.npy'), allow_pickle=True)
		Fneu = np.load(folder.joinpath('Fneu.npy'), allow_pickle=True)
		iscell = np.load(folder.joinpath('iscell.npy'), allow_pickle=True)

		return F, Fneu, stat, ops, iscell

	@staticmethod
	def _load_stat_from_folder(folder):
		return np.load(folder.joinpath('stat.npy'), allow_pickle=True)

	@staticmethod
	def _load_ops_from_folder(folder):
		return np.load(folder.joinpath('ops.npy'), allow_pickle=True).item()

	@staticmethod
	def _load_Fneu_from_folder(folder):
		return np.load(folder.joinpath('Fneu.npy'), allow_pickle=True)

	@staticmethod
	def _load_F_from_folder(folder):
		return np.load(folder.joinpath('F.npy'), allow_pickle=True)

	@staticmethod
	def _load_iscell_from_folder(folder):
		return np.load(folder.joinpath('iscell.npy'), allow_pickle=True)


class Suite2pMultiPlaneLoader(Suite2pLoader):
	"""
	For loading data from multiple plane folders (useful for roi pixel locations)

	:param: load_folder: Path loadfolder: parent suite2p directory, which holds all the plane0, plane1,...subdirs.
	:param: str expected_folder_name: pattern expected for plane folders (default="plane*")

	"""
	expected_folder_name = "suite2p*"
	Suite2pMultiPlaneOutputs = namedtuple("Suite2pMultiPlaneOutputs", ["F", "Fneu", "stat", "ops", "iscell"])

	def __init__(self, **kwargs):

		if "load_folder" in kwargs:  # initialize loader w/ path directly to suite2p/combined
			self.load_folder = Path(kwargs['load_folder'])

		if "good_plane_numbers" in kwargs:
			self.good_plane_numbers = kwargs['good_plane_numbers']

	@property
	def load_folder(self):
		""" Main suite2p folder """
		return self._load_folder

	@load_folder.setter
	def load_folder(self, value):
		"""
		Sets the top level suite2p folder, and automatically updates good_plane_numbers to include every
		subfolder starting with 'plane'
		"""
		self._load_folder = value
		self._all_plane_folders = all_plane_folders(value)
		self.good_plane_numbers = plane_number(self._all_plane_folders)

	@property
	def good_plane_numbers(self) -> list[int]:
		"""List of ints, indicating which planes to load from."""
		return self._good_plane_numbers

	@good_plane_numbers.setter
	def good_plane_numbers(self, value: list[int]):
		"""Sets good_plane_numbers, and automatically updates plane_folders"""
		self._good_plane_numbers = value
		self._plane_folders = list(map(lambda x: self.load_folder.joinpath(f"plane{x}"),
		                               self.good_plane_numbers))

	@property
	def plane_folders(self) -> list[Path]:
		"""Returns list of plane folders corresponding to self.good_plane_numbers"""
		return self._plane_folders

	def load_good_planes(self, **kwargs):
		"""
		loads from selected plane folders

		:param good_plane_numbers : list of ints corresponding to plane folders.
		:return
			- list[Suite2pOutputs(...), Suite2pOutputs(...), ]
				List of namedtuple Suite2pOutputs (1 per plane)
		"""

		if "good_plane_numbers" in kwargs:
			self.good_plane_numbers = kwargs['good_plane_numbers']

		multiplane_data = []
		for f in self.plane_folders:
			F, Fneu, stat, ops, iscell = Suite2pLoader.load_from_folder(f)
			multiplane_data.append(self.Suite2pOutputs(F, Fneu, stat, ops, iscell))

		return multiplane_data

	def load_var_from_folders(self, value):
		"""
			Specify which variable to load from folders

			Example
			-------
			>>> mploader = Suite2pMultiPlaneLoader(load_folder=folder)
			>>> print(mploader.good_plane_numbers)

				Out[]: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]


			>>> stat_list = list(mploader.load_var_from_folders('stat'))
			>>> plane_idx=0; cell_idx=10;
			>>> type(stat_list[plane_idx][cell_idx])

				Out[]: dict

		:param value: variable/datafile to load from suite2p/plane{#}
					  {'F', 'Fneu', 'iscell', 'stat', 'ops'}
		:return:
		"""
		if value == 'F':
			func = Suite2pLoader._load_F_from_folder
		elif value == 'Fneu':
			func = Suite2pLoader._load_Fneu_from_folder
		elif value == 'iscell':
			func = Suite2pLoader._load_iscell_from_folder
		elif value == 'stat':
			"""  Returns a list of np.ndarrays; stat_list[plane][cell_idx] is a dict          """
			func = Suite2pLoader._load_stat_from_folder
		elif value == 'ops':
			func = Suite2pLoader._load_ops_from_folder

		dirlist = map(lambda x: self.load_folder.joinpath(f"plane{x}"), self.good_plane_numbers)

		return list(map(func, dirlist))

	def load(self, **kwargs):
		"""
		Loads F, Fneu, stat, ops, and iscell from every plane folder.

		Set planes in good_plane_numbers.

		:param groupby: {'var', 'plane'}
		:return:
				If groupby=='plane':
					Returns list of namedtuple Suite2pOutputs(F, Fneu, stat, ops, iscell),
					with length = # of self.good_plane_numbers

				Elif groupby='var':
					returns list of lists, of the form

					[[F_0, F_1, ...], [Fneu_0, Fneu_1,...], [stat_0, stat_1,...],
						[ops_0, ops_1,...], [iscell_0, iscell_1, ...]]
		"""
		if 'groupby' not in kwargs:
			groupby = 'var'

		if groupby == 'plane':
			multiplane_data = self.load_good_planes(**kwargs)
		elif groupby == 'var':
			varlist = ['F', 'Fneu', 'stat', 'ops', 'iscell']
			multiplane_data = [self.load_var_from_folders(item) for item in varlist]
			multiplane_data = Suite2pMultiPlaneLoader.Suite2pMultiPlaneOutputs(*multiplane_data)

		return multiplane_data


class Suite2pCombinedLoader(Suite2pLoader):
	""" Data loader for suite2p outputs from <movie>/**/suite2p/combined """

	expected_folder_name = 'combined'

	def __init__(self, **kwargs):

		if "check_folder_name" in kwargs:
			self.check_folder_name = kwargs['check_folder_name']
		else:
			self.check_folder_name = True

		if "load_folder" in kwargs:  # initialize loader w/ path directly to suite2p/combined
			self.load_folder = Path(kwargs['load_folder'])
		elif "statfile" in kwargs:  # initialize loader w/ filepath to something in the suite2p/combined/... folder
			self.load_folder = Suite2pLoader.statfile_parent(kwargs['statfile'])

	@property
	def load_folder(self):
		""" Folder to load .npy files from """
		return self._load_folder

	@load_folder.setter
	def load_folder(self, value):
		"""
		:param Union[str, Path] value: directory path
		"""
		if self.check_folder_name:
			if value.name == self.expected_folder_name:
				self._load_folder = value
			else:
				raise ValueError(f"name of Suite2pCombinedPlaneLoader.plane0_dir must be {self.expected_folder_name} ")
		else:  # any folder name is fine
			self._load_folder = value

	def load(self):  # returns data from suite2p outputs, unprocessed, in
		""" Loads basic suite2p output files.

		:return
			namedtuple: Suite2pOutputs(F, Fneu, stat, ops, iscell) if all files exist

		:example:
		>>>     s2p_outputs = combined_plane_loader.load()
		>>>     Fc = s2p_outputs.F - 0.7 * s2p_outputs.Fneu

		"""

		F, Fneu, stat, ops, iscell = Suite2pLoader.load_from_folder(self.load_folder)
		return self.Suite2pOutputs(F, Fneu, stat, ops, iscell)


class Suite2pPlane0Loader(Suite2pLoader):  # load single plane suite2p data
	""" Single plane data loader """

	expected_folder_name = 'plane0'

	def __init__(self, **kwargs):
		if "check_folder_name" in kwargs:
			self.check_folder_name = kwargs['check_folder_name']
		else:
			self.check_folder_name = True

		if "load_folder" in kwargs:
			self.load_folder = Path(kwargs['load_folder'])
		elif "statfile" in kwargs:
			self.load_folder = Suite2pLoader.statfile_parent(kwargs['statfile'])

	@property
	def load_folder(self):
		""" Folder to load .npy files from """
		return self._load_folder

	@load_folder.setter
	def load_folder(self, value):
		""" Set load_folder, w/ folder name checking """
		if self.check_folder_name:
			if value.name == self.expected_folder_name:
				self._load_folder = value
			else:
				raise ValueError(f"name of Suite2pPlane0Loader.load_folder must be {self.expected_folder_name} ")
		else:  # any folder name is fine
			self._load_folder = value

	def load(self):
		""" Loads basic suite2p output files."""
		F, Fneu, stat, ops, iscell = Suite2pLoader.load_from_folder(self.load_folder)
		return self.Suite2pOutputs(F, Fneu, stat, ops, iscell)


# %%
def main():
	# Testing Suite2pCombinedPlaneLoader() class
	statfile = r".\data\processed_data\2021-08-30\4\movie_003\suite2p\combined\stat.npy"
	combined_plane_loader = Suite2pCombinedLoader(statfile=Path(statfile))
	soutputs = combined_plane_loader.load()
	print(f"\ncombined_plane_loader: ")
	print("---------------------")
	pp.pprint(vars(combined_plane_loader))
	print(f"s2p_outputs = {soutputs.__doc__}")

	combdir = r".\data\processed_data\2021-08-30\4\movie_003\suite2p\combined"
	combined_plane_loader = Suite2pCombinedLoader(load_folder=Path(combdir))
	soutputs = combined_plane_loader.load()
	print(f"\ncombined_plane_loader: ")
	print("---------------------")
	pp.pprint(vars(combined_plane_loader))
	print(f"s2p_outputs = {soutputs.__doc__}")

	# Testing Suite2pCombinedPlaneLoader() class

	combdir = r".\data\processed_data\2021-06-16\1\movie_002\suite2p\plane0"
	plane0_loader = Suite2pPlane0Loader(load_folder=Path(combdir))
	s2p_outputs = plane0_loader.load()
	print(f"\nplane0_loader: ")
	print("-----------------")
	pp.pprint(vars(plane0_loader))
	print(f"plane0_loader = {plane0_loader.__doc__}")

	statfile = r".\data\processed_data\2021-08-30\2\movie_001\1\suite2p\plane0\stat.npy"
	plane0_loader = Suite2pPlane0Loader(statfile=Path(statfile))
	s2p_outputs = plane0_loader.load()
	print(f"\nplane0_loader: ")
	print("-----------------")
	pp.pprint(vars(plane0_loader))
	print(f"plane0_loader = {plane0_loader.__doc__}")

	# load multiplane data
	suite2p_folder = Path(r".\data\processed_data\2021-08-30\4\movie_003\suite2p")
	multiplane_loader = Suite2pMultiPlaneLoader(load_folder=suite2p_folder)
	multiplane_loader.good_plane_numbers = list(range(0, 16, 4))
	s2p_outputs = multiplane_loader.load()

	blk_meta_loader = BlockMetaDataLoader.init_from_statfile(stat_file_1)
	myaml = blk_meta_loader.load_metadata_yaml()
	session_info = blk_meta_loader.load_session_info()
	sync_meta = blk_meta_loader.load_sync_meta()
	xml_meta = blk_meta_loader.load_xml_meta()
	# %%
	stat_file_2 = Path(
		r"G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\data\processed_data\2021-08-30\4\movie_001\suite2p\combined\stat.npy")
	mov_meta_loader = MovieMetaDataLoader.init_from_statfile(stat_file_2)
	myaml = mov_meta_loader.load_metadata_yaml()
	session_info = mov_meta_loader.load_session_info()
	sync_meta = mov_meta_loader.load_sync_meta()
	xml_meta = mov_meta_loader.load_xml_meta()


if __name__ == '__main__':
	main()
