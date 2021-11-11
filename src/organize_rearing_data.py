"""
Write metadata files for rearing datasets in an organized, numbered form in a different folder.
./mb-odor-rearing/data/processed/<fly_uid>/
"""
import pandas as pd
import numpy as np
from pathlib import Path
import json
import collections

import shutil

import rearing_datasets
from ryim import dirutils
import step00_process_suite2p

NEAT_DIR = Path("../mb-odor-rearing/data/processed").resolve()

db = pd.read_csv("src/db.csv")
db.set_index(['date', 'fly_num'], inplace=True)

pfo_ds = rearing_datasets.get_pfo_datasets()
odor_ds = rearing_datasets.get_odor_datasets()
rearing_ds = pfo_ds + odor_ds

pfo_stat_files = rearing_datasets.get_pfo_statfiles()
odor_stat_files = rearing_datasets.get_pfo_statfiles()
stat_files = pfo_stat_files + odor_stat_files


#%%
Path.cwd()

fly_manifest = []
for ds in rearing_ds:
	flydate = ds.flydate
	flynum = int(ds.flynum)
	rear_type = ds.rear_type
	fly_uid = int(db.loc[(flydate, flynum)]['fly_uid'])

	flydict = collections.OrderedDict(fly_uid=fly_uid,
	                                  flydate=flydate,
	                                  flynum=flynum,
	                                  rear_type=rear_type)
	fly_manifest.append(flydict)

	flydir = NEAT_DIR.joinpath(str(fly_uid))
	flydir.mkdir(parents=True, exist_ok=True)
	with open(flydir.joinpath('fly_info.json'), 'w') as out_file:
		json.dump(flydict, out_file, indent=4)

	# movie list
	mov_stat_files = rearing_datasets.path_to_statfiles(ds)
	mov_list = ds.mov_list

	for i, (statfile, mov) in enumerate(zip(mov_stat_files, mov_list)):what
		ainfo = dirutils.parse_dataset_info(statfile)

		movdir = flydir.joinpath(str(i))
		movdir.mkdir(parents=True, exist_ok=True)

		nt_meta_data, nt_s2p_data, title_str = step00_process_suite2p.load_all(statfile)

		nt_meta_data.df_stimulus.to_csv(movdir.joinpath('df_stimulus.csv'), index=False)
		nt_meta_data.df_frametimeinfo.to_csv(movdir.joinpath('df_frametimeinfo.csv'), index=False)
		np.save(movdir.joinpath('timestamps.npy'), nt_meta_data.timestamps)

		xml_meta = nt_meta_data.xml_meta
		with open(movdir.joinpath('thorimage_meta.json'), 'w') as out_file:
			json.dump(xml_meta, out_file, indent=4)

		sync_meta = nt_meta_data.sync_meta
		for key, value in sync_meta.items():
			if isinstance(value, np.ndarray):
				sync_meta[key] = value.tolist()
		with open(movdir.joinpath('sync_meta.json'), 'w') as out_file:
			json.dump(sync_meta, out_file, indent=4)

		session_info = nt_meta_data.session_info
		session_info['rearing'] = rear_type
		mov_info = collections.OrderedDict(movnum=i,
		                                   name=mov.name,
		                                   chunk=mov.chunk,
		                                   conc=mov.concentration,
		                                   planes=mov.planes,
		                                   frames=mov.frames,
		                                   nacq=mov.nacq,
		                                   quality=mov.quality,
		                                   microscope=nt_meta_data.session_info['microscope'],
		                                   thorimagename=nt_meta_data.session_info['thorimagename'],
		                                   thorsyncname=nt_meta_data.session_info['thorsyncname'],
		                                   )
		with open(movdir.joinpath('mov_info.json'), 'w') as out_file:
			json.dump(mov_info, out_file, indent=4)

		# copy trial peak outputs into organized directory
		source_dir = statfile.with_name('step00_process_suite2p')
		destination_dir = movdir.joinpath("step00_process_suite2p")
		shutil.copytree(source_dir, destination_dir, dirs_exist_ok=True)

		source_dir = statfile.with_name('step01_extract_peaks')
		destination_dir = movdir.joinpath("step01_extract_peaks")
		shutil.copytree(source_dir, destination_dir, dirs_exist_ok=True)

		source_dir = statfile.with_name('step02_extract_trial_tensors')
		destination_dir = movdir.joinpath("step02_extract_trial_tensors")
		shutil.copytree(source_dir, destination_dir, dirs_exist_ok=True)

#%%

# save full fly manifest
with open(NEAT_DIR.joinpath('fly_manifest.json'), 'w') as out_file:
	json.dump(fly_manifest, out_file, indent=4)

# only odor-reared flies
odor_reared_manifest = [item for item in fly_manifest if item['rear_type'] == '1-6ol+EP']
odor_reared_manifest = sorted(odor_reared_manifest, key=lambda x: x['fly_uid'])

with open(NEAT_DIR.joinpath('odor_reared_manifest.json'), 'w') as out_file:
	json.dump(odor_reared_manifest, out_file, indent=4)

# only pfo-reared flies
pfo_reared_manifest = [item for item in fly_manifest if item['rear_type'] == 'pfo']
pfo_reared_manifest = sorted(pfo_reared_manifest, key=lambda x: x['fly_uid'])

with open(NEAT_DIR.joinpath('pfo_reared_manifest.json'), 'w') as out_file:
	json.dump(pfo_reared_manifest, out_file, indent=4)











