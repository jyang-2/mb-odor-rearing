from ryim import thorimagels
from pathlib import *

# %% STEP 1: MAKE METADATA.YAML
"""
- Check metadata.yaml in SublimeText
- Run make_excel_metadata.py to create flyinfo.xlsx (in every fly folder) and metadata.xlsx (in every movie folder)
- Output files: 
    - flyinfo.xlsx in each fly folder
    - metadata.xlsx in each movie folder
     
"""


# %% STEP 2: Create yaml and json info files
"""
    Use thorimagels 
"""
raw_data_folder = Path(r"G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\data\2p\2021-06-16") # folder in mb_odor_rearing/data/2p
proc_data_folder = Path(r"G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\data\processed_data\2021-06-16")
xml_files = sorted(list(raw_data_folder.glob("**/Experiment.xml")))

# first, run in data/2p folders
#%% STEP 3: PROCESS EXPERIMENT.XML FILES --> XML_META.YAML and XML_META.JSON

print("\nParsing Experiment.xml files in data/2p")
for f in xml_files:
    print(f)
    thorimagels.move_to_processed_data(f)  # copy Experiment.xml to data/processed folder
    thorimagels.xml_to_yaml(f)
    thorimagels.xml_to_json(f)
    print('----------DONE-----------')

print("\nMoving xml_meta.yaml files to data/processed_data")
for f in sorted(list(raw_data_folder.glob("**/xml_meta.yaml"))):
    print(f"\t{f}")
    thorimagels.move_to_processed_data(f)

print("\nMoving xml_meta.json files to data/processed_data")
for f in sorted(list(raw_data_folder.glob("**/xml_meta.json"))):
    print(f"\t{f}")
    thorimagels.move_to_processed_data(f)


#%% STEP 3: MATLAB - run_syncdata_to_json.m
# /src/matlab/huginn - run Untitled.m (or run_syncdata_to_json.m if the name was changed)

print(f"\nCopying sync_data.json files from data/processed_data to data/2p (for movie conversion info)")
print(f"----------------------------------------------------------------------------------------------")
json_files = sorted(list(proc_data_folder.glob("**/sync_meta.json")))
for f in json_files:
    print(f"\t{f}")
    thorimagels.move_to_2p_data(f)
    print("\t\tcopied")

#%% STEP 4 : CONVERT .RAW MOVIES TO TIFF STACKS

print(f"\nConverting .raw movies into tiff stacks in data/2p")
print(f"----------------------------------------------------")

raw_files = sorted(list(raw_data_folder.glob("**/Image_*.raw")))
for f in raw_files:
    print(f"\n\t{f}")
    thorimagels.dask_raw2tiffs(f)
    print(f"\t\tdone")

print(f"\nMoving tiff files from data/2p --> data/processed data")
print(f"----------------------------------------------------")

raw_tiffs = sorted(list(raw_data_folder.glob("**/raw_tifs/stk*.tif")))
for f in raw_tiffs:
    print(f"\n\t{f}")
    thorimagels.move_to_processed_data(f)
    print(f"\t\tdone")

print(f"\nMoving tiff files in raw_tifs folder up 1 directory")
print(f"----------------------------------------------------")

moved_tiffs = sorted(list(proc_data_folder.glob("**/raw_tifs/stk*.tif")))
for f in moved_tiffs:
    print(f"\n\t{f}")
    thorimagels.up_one_dir(f)
    print(f"\t\tdone")

#%% RUN SUITE 2P, THEN ANALYZE RESULTS

#%%
rearing
