import pandas as pd
import numpy as np
from natsort import natsorted


def get_chunksubdirs(folder):
    """

    Parameters
    ----------
    folder : directory (containing full-movie associated files like 'metadata.xlsx') to search for chunks (split movies)

    Returns
    -------
    chunksubdirs : List(str) of subdirectories w/ integer-only names

    """
    subdirs = [item.name for item in folder.iterdir() if item.is_dir()]
    chunksubdirs = natsorted(filter(str.isnumeric, subdirs))
    return chunksubdirs


def split_timestamps(fname):
    subdirs = get_chunksubdirs(fname.parent)
    n = len(subdirs)
    if n > 0:
        save_files = map(lambda x: fname.parent.joinpath(x, 'timestamps.npy'), subdirs)
        save_files = list(save_files)

        timestamps = np.load(fname)
        ts_list = np.split(timestamps, n)
        for f, ts in zip(save_files, ts_list):
            np.save(f, ts)
        return save_files
    else:
        return None


def split_metadata_xlsx(fname):
    subdirs = get_chunksubdirs(fname.parent)
    n = len(subdirs)

    if n > 0:
        metadata = pd.read_excel(fname, sheet_name='stimulus')
        chunksize = metadata.shape[0] / n
        metadata['chunk'] = np.repeat(range(n), chunksize)

        saved_files = []
        grouped_df = metadata.groupby(by='chunk')
        for key, item in grouped_df:
            savename = fname.parent.joinpath(str(key), 'metadata.csv')
            grouped_df.get_group(key).to_csv(savename, index=False)
            saved_files.append(savename)
        return saved_files
    else:
        return []


def split_frametimeinfo(fname, edit_original=True):
    n = len(get_chunksubdirs(fname.parent))
    frametimeinfo = pd.read_csv(fname)
    chunksize = frametimeinfo.shape[0] / n
    frametimeinfo['chunk'] = np.repeat(range(n), chunksize)

    if edit_original:
        frametimeinfo.to_csv(fname, index=False)

    grouped_df = frametimeinfo.groupby(by='chunk')
    for key, item in grouped_df:
        savename = fname.parent.joinpath(str(key), 'frametimeinfo.csv')
        grouped_df.get_group(key).to_csv(savename, index=False)

    return frametimeinfo

    # %%
    fti_list = np.array_split(frametimeinfo, n)


if __name__ == '__main__':
    if chunkdirs is not None:
        if len(chunkdirs) != n:
            raise ValueError('chunkdirs must be of length n')

        for fti, chunk in zip(fti_list, chunkdirs):
            if isinstance(chunk, str):
                chunk = int(chunk)
            fti['chunk'] = chunk

    frametimeinfo = pd.concat(fti_list)
