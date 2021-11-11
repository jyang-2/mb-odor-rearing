HALLEM.MAT CONTENTS
	receptor
		all: [1x24 string], names of all glomeruli
		phr: [1x4 string], names of pheremone glomeruli
	     nonphr: [1x20 string], names of non-pheremone glomeruli


	chem
	    cat_list: [10x1 string], chemical categories
	      	cmap: [10x3 double], colormap of chem. categories, out of 255
	    cmap_rgb: [10x3 double], colormap of chem. categories, out of 1


	odors
	   odor_list: [111x1 string], names of hallem odors
	      abbrev: [111x1 string], abbreviations for hallem odors

    pairs
        A: [111 x 111 matrix], index into with (i,j) for hpair_ind
        T: [12321 x 3] table, with variables hind1, hind2, hpair_ind