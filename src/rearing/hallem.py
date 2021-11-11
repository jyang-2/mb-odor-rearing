from pathlib import Path
import numpy as np
import pandas as pd


# colormap: xkcd_rgb
# amines: gold
# lactones: indigo
# acids: magenta
# sulfur compounds: black
# terpenes: green/yellow
# aldehydes: gray
# ketones: yellow
# aromatics: light blue
# alcohols: red
# esters: dark green


class HallemSet:
    def __init__(self, start_index=1):
        filepath = Path(__file__).parent
        df = pd.read_csv(filepath.joinpath("hallem-8ea3408c-de8c-42a6-a609-8bf4c3770264.csv"))
        df = df.set_index('odor')
        self.df = df

        # odors = dict()
        # odors['odor_list'] = df.index.tolist()
        # odors['abbrev_list'] = df.abbrev.tolist()
        # odors['hi_list'] = df.hallem_ind.tolist()
        # self.odors = odors
        self.odor_list = df.index.tolist()
        self.abbrev_list = df.abbrev.tolist()
        self.hi_list = df.hallem_ind.tolist()

        # receptor = dict()
        receptors_all = df.columns.tolist()[:24]
        receptors_phr = ["x33b", "x47b", "x65a", "x88a"]
        receptors_nonphr = [item for item in receptors_all if item not in receptors_phr]
        self.receptors_all = receptors_all
        self.receptors_phr = receptors_phr
        self.receptors_nonphr = receptors_nonphr

        self.chemcat_list = ["amines", "lactones", "acids", "sulfur compounds", "terpenes", "aldehydes", "ketones",
                             "aromatics", "alcohols", "esters"]

        self.chemcat_cmap = [[229, 183, 43],
                             [56, 17, 127],
                             [229, 66, 188],
                             [0, 0, 0],
                             [135, 216, 95],
                             [160, 160, 160],
                             [255, 246, 94],
                             [191, 236, 239],
                             [191, 25, 0],
                             [34, 94, 29]]

    def odor2hi(self, odor, usenan=False):
        if usenan:
            valid_hallem_odor = odor in set(self.odor_list)
            if valid_hallem_odor:
                idx = self.odor_list.index(odor)
                hi = self.hi_list[idx]
            else:
                hi = -1
        else:
            idx = self.odor_list.index(odor)
            hi = self.hi_list[idx]
        return hi

    def odor2abbrev(self, odor):
        idx = self.odor_list.index(odor)
        abbrev = self.abbrev_list[idx]
        return abbrev

    def hi2odor(self, hi):
        idx = self.hi_list.index(hi)
        odor = self.odor_list[idx]
        return odor

    def hi2abbrev(self, hi):
        idx = self.hi_list.index(hi)
        abbrev = self.abbrev_list[idx]
        return abbrev

    def hipair2str(self, hipair):
        hi1, hi2 = hipair
        s = set(hipair)
        if len(s) == 1:
            idx = self.hi_list.index(hi1)
            st = self.abbrev_list[idx]
            return st
        elif s == {111, 112}:
            idx = self.hi_list.index(111)
            st = self.abbrev_list[idx]
            return st
        else:
            s.discard(111)
            s.discard(112)
            x = [self.hi2abbrev(item) for item in s]
            st = '+'.join(x)
            return st

    def PNtransform(self, ORNc):
        # parameters for transformation
        rmax = 165
        n = 1.5
        sigma = 12
        m_inp = 10.63

        # zero all negative firing rates
        if isinstance(ORNc, pd.DataFrame):
            ORN = ORNc.to_numpy()
        ORN[ORN < 0] = 0

        # compute the EAG based on raw Hallem ORN values
        EAG = ORN.sum(axis=1)
        EAG = EAG[:, np.newaxis]
        EAG[EAG < 0] = 0

        # compute suppression factor
        s = m_inp * EAG / 190

        # input gain, homogeneous, set at VM7 levels from integration study
        ipn_homo = rmax * (ORN ** n) / ((ORN ** n) + sigma ** n + s ** n)
        if isinstance(ORNc, pd.DataFrame):
            return pd.DataFrame(ipn_homo, index=ORNc.index, columns=ORNc.columns)
        else:
            return ipn_homo


def main():
    print(__file__)


if __name__ == "__main__":
    main()
