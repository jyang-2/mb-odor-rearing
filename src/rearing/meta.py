import pandas as pd
from rearing import hallem
from itertools import product

Hallem = hallem.HallemSet()

def get_pin_odor_conc(s):
    """
    Parses string in odor_panel (from stimulus_info.yaml) to find the odor concentration.

    :param str s: String similar to "<odor> (<log_10_conc>)" to parse for concentration value
    :return: conc - log10 odor concentration
    :rtype: float
    """
    conc = s[s.find("(") + 1:s.find(")")]
    conc = float(conc)
    return conc


def get_pin_odor(s):
    """
    Parses string in odor_panel (from stimulus_info.yaml) to find the odor name.

    :param str s: String similar to "<odor> (<log_10_conc>)" to parse for odor name (full, hallem format)
    :return: odor_name
    :rtype: str
    """
    odor_name = s[0:s.find("(")]
    odor_name = str.strip(odor_name)
    return odor_name


def get_session_info(acquisitions, exp):
    """
    Makes dataframe w/ columns ['date', 'fly', 'microscope', 'thorimagename', 'thorsyncname'] to put in 'sessioninfo'
    sheet of metadata.xlsx. Contains list of thorimage and thorsync names, along with experiment information.

    example:
        expt = dict()
        expt['date'] = '2021-04-23'
        expt['fly'] = 2
        expt['microscope'] = 'galvo-res'

        metadata = io.load_yaml('metadata.yaml')
        acq = metadata['Acquisitions']['fly2']

        df_sessioninfo = meta.get_sessioninfo(acquisitions, expt)

    :param acquisitions: list of dicts describing image acquisitions, acquisitions = metadata['Acquisitions'][flyID]
    :type acquisitions: list
    :param exp:
    :type exp: dict
    :return: df_sessioninfo, DataFrame w/ columns ['date', 'fly', 'microscope', 'thorimagename', 'thorsyncname']
    :rtype: pd.DataFrame
    """
    # make sessioninfo dataframe
    df_session_info = pd.DataFrame(
        [(exp['date'], exp['fly'], exp['microscope'], acq['thorimage'], acq['thorsync'], acq['notes']) for acq in acquisitions],
        columns=['date', 'fly', 'microscope', 'thorimagename', 'thorsyncname', 'notes'])
    return df_session_info


def get_channel_info(channel, df_olf):

    """
    metadata = io.load_yaml('metadata.yaml')
    acquisitions = metadata['Acquisitions']['fly2']

    odor_panel = metadata['Odor delivery']['pin odors'][1]
    df_olf = meta.get_olfactometer_info(odor_panel)

    pin_list = acquisitions[0]['stimulus']['channelA']
    chanA = get_channel_info(pin_list, df_olf)

    :param channel: list of pins used in arduino stimulus delivery, ex. [2,5,6]
    :type channel: list

    :param df_olf: Dataframe of olfactometer information (columns = pin, hi, odor, odor_abbrev, conc)
    :type df_olf: pandas.DataFrame

    :return : df_channel - Dataframe with olfactometer information (pin, hi, odor, odor_abbrev, conc)
                           in stimulus channel order.
    :rtype: pandas.DataFrame
    """
    df_olf = df_olf.set_index('pin')
    df_channel = df_olf.loc[channel]
    df_channel = df_channel.reset_index()
    return df_channel


def get_olfactometer_info(panel, as_dataframe=True):
    """
    example:
        metadata = io.load_yaml('metadata.yaml')
        odor_panel = metadata['Odor delivery']['pin odors'][1]
        df_olf = meta.get_olfactometer_info(odor_panel)

    :param panel: dictionary of pin-odor mappings defined in metadata.yaml
    :type panel: dict
    :param as_dataframe: whether to return olfactometer variables as a dataframe or separate vars
    :type as_dataframe: boolean
    :return: df_olf, with columns {'pin', 'hi', 'odor', 'conc'}, or odor_panel_trim, odor_list, hi_list, conc_list, pin_list
    :rtype: pd.DataFrame


    """
    odor_panel_trim = {x: y for x, y in panel.items() if y is not None}
    odor_list = [get_pin_odor(s).lower() for s in odor_panel_trim.values()]
    hi_list = [Hallem.odor2hi(item) for item in odor_list]
    conc_list = [get_pin_odor_conc(s) for s in odor_panel_trim.values()]
    pin_list = [k for k in odor_panel_trim.keys()]

    if as_dataframe:
        df_olf = pd.DataFrame({'pin': pin_list,
                               'hi': hi_list,
                               'odor': odor_list,
                               'conc': conc_list})
        return df_olf
    else:
        return odor_panel_trim, odor_list, hi_list, conc_list, pin_list


def get_stimtype(hi_list):
    hipairs = list(product(hi_list, hi_list))
    hi1, hi2 = zip(*hipairs)
    stimname = [Hallem.hipair2str(item) for item in hipairs]
    df_stimtype = pd.DataFrame({'hi1': hi1, 'hi2': hi2, 'stimname': stimname})
    return df_stimtype
