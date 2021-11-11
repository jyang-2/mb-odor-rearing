import numpy as np
from scipy import spatial, signal, stats, ndimage, optimize
from sklearn import preprocessing


def quantile_match(moving, ref,  nbins, type='non-negative'):
