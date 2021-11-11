import numpy as np
from pathlib import Path
from tempfile import TemporaryDirectory
import yaml
import json
import tifffile

from suite2p.registration import register
import suite2p

refImg = register.pick_initial_reference(ops)

ops = suite2p.run_s2p.default_ops()

ops_file = Path.cwd().joinpath('data', 'processed_data', '2021-08-30', '2', 'movie_002', 'ops.npy')
ops.update(np.load(ops_file, allow_pickle=True).item())
ops.update({'data_path': ['ops.parent'],
            'save_path0': TemporaryDirectory().name,
            'tiff_list': ['stk_0000.tif']
            })



ops.parent.joinpath()

yaml_file = raw_file.parent.joinpath("xml_meta.yaml")
with open(yaml_file, 'r') as f:
    meta = yaml.safe_load(f)
