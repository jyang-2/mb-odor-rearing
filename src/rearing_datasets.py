"""
Information about available datasets, parameters used in rest of project, and useful functions.

ex:
---

    pfo_datasets = rearing_datasets.get_pfo_datasets()
    odor_stat_files = rearing_datasets.get_odor_statfiles()
"""

from pathlib import *
from ryim import dirutils, analysis
from collections import namedtuple
import itertools

odor_order = ['paraffin', 'pentanoic acid', 'benzaldehyde', '1-pentanol', '1-hexanol', 'ethyl propionate']


clr_dict = {
    'pfo': 'gray',
    '1-5ol': 'tab:blue',
    '1-6ol': 'tab:green',
    'EP': 'tab:red',
    'VA': 'tab:brown',
    '1-6ol+EP': 'tab:purple',
}

DataSet = namedtuple("DataSet", ["flydate", "flynum", "rear_type", "mov_list"])

MovieSet = namedtuple("MovieSet", ["name", "chunk", "concentration", "planes", "frames", "nacq", "quality"])
# if chunk = -1, it is not a chunked movie
# concentration > 0, then the movie contains a mix of different concentrations.

pfo_datasets = [DataSet('2021-08-30', 2, 'pfo', [MovieSet("movie_001", 0, -5, 16, 2907, 3, 'good'),
                                                 MovieSet("movie_001", 1, -4, 16, 2907, 3, "good"),
                                                 MovieSet("movie_001", 2, -3, 16, 2907, 3, "good"),
                                                 MovieSet("movie_002", -1, -3, 16, 2907, 3, "good")
                                                ]
                        ),
                DataSet('2021-08-24', 1, 'pfo', [MovieSet("movie_002", -1, -3, 16, 2907, 3, "ok-good"),
                                                 MovieSet("movie_003", -1, -5, 16, 2907, 3, "ok-good"),
                                                 MovieSet("movie_004", -1, -4, 16, 2907, 3, "ok-good"),
                                                 MovieSet("movie_005", -1, -3, 16, 2907, 3, "ok-good"),
                                                 MovieSet("movie_006", -1, 1, 16, 2907, 3, "REVISIT")  # bad
                                                 ]
                        ),
                DataSet('2021-06-26', 2, 'pfo',
                        [MovieSet("movie_001", -1, -3, 1, 5796, 3, 'good'),  # win=45, prctile=45
                         MovieSet('movie_002', -1, -5, 1, 5796, 3, 'CORRMAT LOOKS BAD-FIX '),
                         # paraffin response too strong?
                         MovieSet('movie_003', -1, -4, 1, 5796, 3, 'might be okay?'),  # win60, prctile20, peakwin15
                         ]
                        ),
                # young fly
                DataSet('2021-07-27', 1, 'pfo', [MovieSet('movie_001', -1, -3, 16, 2686, 3, "good"),
                                                 MovieSet('movie_002', -1, -5, 16, 2269, 3,
                                                          "corrmat iffy, worth a try"),
                                                 MovieSet('movie_003', -1, -4, 16, 2268, 3, "good"),
                                                 MovieSet("movie_004", -1, -3, 16, 2268, 3, "might be dying, REVIST")
                                                 # looks okay, actually? check traces
                                                 ]
                        ),
                DataSet('2021-06-16', 1, 'pfo', [MovieSet("movie_002", -1, -4, 1, 3440, 3, "some motion/drift"),
                                                 MovieSet('movie_003', -1, -3, 1, 3441, 3, "some motion/drift"),
                                                 MovieSet("movie_004", -1, -5, 1, 3440, 3, "might be dying, REVIST")
                                                 ]
                        )
                ]
#  consider redoing traces
odor_datasets = [DataSet('2021-08-30', 1, '1-6ol+EP', [MovieSet("movie_002", 0, -5, 16, 2907, 3, "good"),
                                                       MovieSet("movie_002", 1, -4, 16, 2907, 3, "good"),
                                                       MovieSet("movie_002", 2, -4, 16, 2907, 3, "good")]
                         ),

                 DataSet('2021-08-30', 4, '1-6ol+EP', [MovieSet("movie_001", -1, -5, 16, 2598, 3, "good"),
                                                       MovieSet("movie_002", -1, -4, 16, 2598, 3, "good"),
                                                       MovieSet("movie_003", -1, -4, 16, 2598, 3, "good")
                                                       ]
                         ),
                 # note - this fly doesn't have VA, has benzaldehyde instead (not at all concentrations)
                 # has sparsity analyses done
                 DataSet('2021-06-26', 1, '1-6ol+EP', [
                     MovieSet("movie_002", -1, -5, 1, 3864, 3, "iffy, trying anyway"),
                     MovieSet("movie_005", -1, -3, 1, 2898, 3, "great"),
                     MovieSet("movie_006", -1, -3, 1, 3864, 3, "great")
                 ]
                         ),
                 DataSet('2021-06-14', 1, '1-6ol+EP',
                         [MovieSet("movie_001", -1, -3, 1, 3449, 3, "good"),  # from the galvo-galvo
                          MovieSet('movie_002', -1, -4, 1, 3447, 3, "movie jumpy, corrmat looks good"),
                          MovieSet('movie_003', -1, -5, 1, 3450, 3, "good"),
                          MovieSet('movie_004', -1, -3, 1, 3449, 3, 'good')])
                 ]


def get_pfo_datasets():
    """
    Returns list of DataSet namedtuples for paraffin-reared flies.

    :return: pfo_datasets
    :rtype: list[Path]

    ex:
    ---
        pfo_ds = rearing_datasets.get_pfo_datasets()


    """
    return pfo_datasets


def get_odor_datasets():
    """
    Returns list of DataSet namedtuples for odor-reared flies.

    :return: odor_datasets
    :rtype: list[Path]

    ex:
    ---
        odor_ds = rearing_datasets.get_odor_datasets()

    """
    return odor_datasets


def path_to_statfiles(ds):
    """
    Accepts a DataSet namedtuple, and returns the expected stat.npy filepaths to ds.mov_list

    :param ds: DataSet(flydate, flynum, rear_type, mov_list)'

    ex:
    ---
        pfo_ds = rearing_datasets.get_pfo_datasets()
        stat_files = path_to_statfiles(pfo_ds[0])

    """

    FLY_DIR = Path("./data/processed_data").joinpath(ds.flydate, str(ds.flynum))

    stat_file_list = [None] * len(ds.mov_list)

    for i, m in enumerate(ds.mov_list):
        if m.chunk < 0:  # not a chunked movie
            S2P_DIR = FLY_DIR.joinpath(m.name, 'suite2p')
        elif m.chunk >= 0:
            S2P_DIR = FLY_DIR.joinpath(m.name, str(m.chunk), 'suite2p')

        if m.planes > 1:
            stat_file = S2P_DIR.joinpath("combined", "stat.npy")
        elif m.planes == 1:
            stat_file = S2P_DIR.joinpath("plane0", "stat.npy")
        stat_file_list[i] = stat_file
    return stat_file_list


def get_pfo_statfiles(list_type='flat'):
    """
    Get list of all filepaths to stat.npy files, for paraffin-reared flies datasets listed in rearing_datasets.py

    :param list_type: item in {'flat', 'nested'}
                        - 'flat' : returns paths to stat.npy in flat list
                        - 'nested' : returns list of lists, where statfiles for MovieSet are grouped by DataSet

    :type list_type: str

    ex:
    ---
        pfo_stat_files = get_pfo_statfiles()
    """

    pfo_stat_files_ = [path_to_statfiles(item) for item in pfo_datasets]
    if type == 'nested':
        return pfo_stat_files_
    else:
        return [item for sublist in pfo_stat_files_ for item in sublist]


def get_odor_statfiles(list_type='flat'):
    """
    Get list of all filepaths to stat.npy files, for odor-reared fly datasets listed in rearing_datasets.py

    :param list_type: item in {'flat', 'nested'}
                        - 'flat' : returns paths to stat.npy in flat list
                        - 'nested' : returns list of lists, where statfiles for MovieSet are grouped by DataSet

    :type str

    ex:
    ---
        odor_stat_files = get_pfo_statfiles()
    """

    """"""

    odor_stat_files_ = [path_to_statfiles(item) for item in odor_datasets]
    if type == 'nested':
        return odor_stat_files_
    else:
        return [item for sublist in odor_stat_files_ for item in sublist]


def check_stat_file(statfile):
    """
    Check that dirutils.parse_datast_info(statfile) works, and that metadata can be loaded
     w/ analysis.BlockMetaDataLoader or analysis.MovieMetaDataLoader.

    :param statfile: path to a suite2p/**/stat.npy file.
    :type statfile: Path | str
    :return:
    True if filestructure for statfile is valid w/ loadable metadata
    False if not
    :rtype: bool

    """
    try:
        ainfo = dirutils.parse_dataset_info(statfile)

        # loading metadata
        if ainfo.is_block:
            meta_data_loader = analysis.BlockMetaDataLoader.init_from_statfile(statfile)
        else:
            meta_data_loader = analysis.MovieMetaDataLoader.init_from_statfile(statfile)
        nt_meta_data = meta_data_loader.load()
        return True
    except Exception as error:
        print('Dataset info for statfile cannot be parsed, or metadata cannot be loaded.')
        print(error)
        return False


def main():
    pfo_stat_files_ = get_pfo_statfiles()
    odor_stat_files_ = get_odor_statfiles()
    return pfo_stat_files_, odor_stat_files_


if __name__ == '__main__':
    pfo_stat_files, odor_stat_files = main()

#
#     failed_save = []
#     for i, f in enumerate(odor_stat_files):
#         ops = np.load(f.with_name('ops.npy').absolute(), allow_pickle=True).item()
#         if Path(ops['save_path']).parts[0] == 'local':
#             failed_save.append(f)
#             print(f'FAILED SAVE : {f}')
#         else:
#             save_folder = fr"{dirutils.parse_parent_suite2p(f)}"
#             suite2p.io.save_nwb(save_folder)
#
#
#     nwbfiles = sorted(list(Path.cwd().glob("**/ophys.nwb")))
#
#     with NWBHDF5IO(str(nwbfiles[-1]), 'r') as io:
#         nwbfile_in = io.read()
#         # F = nwbfile_in.processing['ophys']['Fluorescence'].get_roi_response_series('Fluorescence').data
#         F = nwbfile_in.processing['ophys'].get_data_interface('Fluorescence')
#
#     nwbfile_in.objects[F.object_id]
#
#     rois = nwbfile_in.processing['ophys']['ImageSegmentation']['PlaneSegmentation']['voxel_mask']
#     ps = rois.parent
#
#
# #%%
#     io = NWBHDF5IO(str(nwbfiles[-1]), 'r')
#     nwbfile = io.read()

# pfo_stat_files = [path_to_statfiles(item) for item in pfo_datasets]
# odor_stat_files = [path_to_statfiles(item) for item in odor_datasets]
# stat_file_list = [*pfo_stat_files, *odor_stat_files]
# flat_file_list = list(np.concatenate(stat_file_list).flat)
#
# # flatten odor list
# odor_stat_files = [item for sublist in odor_stat_files for item in sublist]
# pfo_stat_files = [item for sublist in pfo_stat_files for item in sublist]
# np.save('odor_stat_files.npy', odor_stat_files)
# np.save('pfo_stat_files.npy', pfo_stat_files)

# for i, f in enumerate(flat_file_list[4:]):
#     print(f"{i}...f{f}")
#     print('==========================================')
#     check_stat_file(f)
#     #run_correlation_analysis.main(f)
#     run_processed_rasterplot.main(f)
#
# # for file_list in stat_file_list:
# #     print(f"\n\n =====FLY======")
# #     for f in file_list[-1]:
# #         check_stat_file(f)
# #         run_correlation_analysis.main(f)
# #         print("done")

#
# def save_dataset_info_to_npz():
#     """
#     Save experiment listings and filepath info to .npz file, for easy looping through different experiments.
#     :return:
#         True if successful
#     """
#     np.savez('stat_files.npz', stat_file_list =stat_file_list, flat_file_list=flat_file_list,
#              pfo_datasets=pfo_datasets, odor_datasets=odor_datasets)
#     print('stat_files.npz was saved successfully.')
#     return True
# %%
# print("------PFO CHECK-------")
# for file_list in pfo_stat_files:
#     print(f"\n\n =====dataset======")
#
#     for f in file_list:
#         check_stat_file(f)
# #%%
# print("------ODOR-REARED CHECK-------")
# for file_list in odor_stat_files:
#     print(f"\n\n =====dataset======")
#
#     for f in file_list:
#         check_stat_file(f)
# # %%
#
#
# def main():
#     """
#
#     Returns
#     -------
#      db : pd.Dataframe, from src/db.csv
#     dbm : pd.Dataframe, exploded movie column from db
#
#     """
#     db = pd.read_csv('src/db.csv')
#
#     # database, but with exploded movie column
#     dbm = copy.deepcopy(db)
#     dbm = dbm.rename(columns={'movies': 'mov'})
#     dbm['mov'] = dbm['mov'].apply(lambda x: x.split(","))
#     dbm = dbm.explode('mov', ignore_index=True)
#     dbm['mov'] = dbm['mov'].apply(lambda x: x.strip())
#     return db, dbm
#
#
# # %%
#
# db, dbm = main()
#
# # dbm['n_chunks'] = np.nan
# # dbm['chunks'] = np.nan
# # dbm['suite2p_folders'] = np.nan
#
# # loop through each movie and find suite2p folders
# list_suite2p_folders = [None] * dbm.shape[0]
#
# for row in dbm.itertuples():
#     mov_dir = Path('.').joinpath('data', 'processed_data', row.date, str(row.fly_num), row.mov)
#     suite2p_folders = list(mov_dir.rglob("**/suite2p"))
#     print(suite2p_folders)
#     print(f"\n{row.Index}...flydate={row.date}, flynum={row.fly_num}, mov={row.mov}")
#
#     list_suite2p_folders[row.Index] = suite2p_folders
#
# dbm['suite2p_folder'] = list_suite2p_folders
# # %%
# # explode suite2p folders
# dbm_chunks = dbm.explode('suite2p_folder', ignore_index=True)
# dbm_chunks['has_suite2p'] = ~dbm_chunks['suite2p_folder'].isnull()
#
# # loop through suite2p folders, and determine if it's in a chunk or not.
# # If it is a chunk, parse which number.
#
# dbm_chunks['is_chunk'] = False
# dbm_chunks['chunk_idx'] = np.nan
#
# for row in dbm_chunks.itertuples():
#     mov_dir = Path('.').joinpath('data', 'processed_data', row.date, str(row.fly_num), row.mov)
#     print(f"\n{row.Index}...flydate={row.date}, flynum={row.fly_num}, mov={row.mov}")
#     print(row.has_suite2p, type(row.suite2p_folder))
#
#     if row.has_suite2p:
#         relpath = row.suite2p_folder.relative_to(mov_dir)
#         relpath_parts = relpath.parts
#         print(f"relpath parts: {relpath.parts}")
#
#         if len(relpath_parts) > 1:
#             dbm_chunks.at[row.Index, 'is_chunk'] = True
#             dbm_chunks.at[row.Index, 'chunk_idx'] = int(relpath_parts[0])
# %%

# if ~np.isnan(row.suite2p_folder):
#    if dirutils.in_block_dir(row.suite2p_folder):
#        dbm_chunks.at[row.Index, 'is_chunk'] = True
# if dirutils.in_block_dir(row.suite2p_folder) & (~np.isnan(row.suite2p_folder)):
#    dbm_chunks.at[row.Index, 'is_chunk'] = True

# if len(suite2p_folders)>0:
#     chunk_folders = [dirutils.in_block_dir(item) for item in suite2p_folders]
#     dbm['n_chunks'] = np.sum(chunk_folders)
#     if chunk_folders >
#
# for f in suite2p_folders:
#     planes = analysis.all_plane_folders(f)
#     n_planes = len(planes)
#     plane_numbers = analysis.plane_number(planes)
#     print(f"\t{f}: {n_planes} planes")
#
#     if n_planes == 1:
#         dbm.at[row.Index, 'planes'] = n_planes
#             stat_file = mov_dir.joinpath('suite2p')
#
#     elif n_planes > 1:
#         dbm.at[row.Index, 'planes'] = n_planes

# movie/suite2p
# movie/#/suite2p
