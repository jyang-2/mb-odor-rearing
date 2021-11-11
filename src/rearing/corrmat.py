import pandas as pd
import numpy as np


"""

- Long to wide format:
  
        dfs.pivot(index=['stim_idx_r', 'stimname_r'],
                  columns=['stim_idx_c', 'stimname_c'],
                  values="value")


- Sort squareform df by index level :

        df = df.sort_index(axis=0, level='stim_idx')


- Sort squareform df by column level :

        df = df.sort_index(axis=1, level='stim_idx')
"""


def stretch(df, mask=None, value_name=None, as_series=True):
    """
    Similar to stretch() function in R package corrr
    :param as_series:
    :type as_series:
    :param value_name:
    :type value_name:
    :param mask:
    :type mask:
    :param df:
    :type df:
    :return:
    :rtype:
    """

    dfs = df.melt(ignore_index=False)

    if mask is not None:
        if mask == "triu":
            keep = np.triu(np.ones(df.shape), k=1)
        elif mask == "tril":
            keep = np.tril(np.ones(df.shape), k=1)
        dfs = dfs[keep.astype('bool').reshape(df.size, order='F')]

    row_names = [item + "_r" for item in list(dfs.index.names)]
    col_names = [item + "_c" if (item != "value") else item for item in dfs.columns.to_list()]

    dfs.index.set_names(row_names, inplace=True)
    dfs.columns = col_names
    dfs.reset_index(inplace=True)

    # rename the value column
    if value_name is not None:
        dfs.rename(columns={"value": value_name}, inplace=True)

    # convert to Series (for combining corr. mat distance metrics)
    if as_series:
        dfs.set_index(['stim_idx_r', 'stim_idx_c', 'stimname_r', 'stimname_c'],
                      inplace=True)

    return dfs


