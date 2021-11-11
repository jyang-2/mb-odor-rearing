"""
Load all the long-form odor-odor correlation tables, and make a summary figure.
"""
from pathlib import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import namedtuple
import pyarrow.feather as feather
import seaborn as sns


# %%
def reflect_tidy_corr(tidy_df):
	"""
    Make tidy odor-odor pair correlation (long format pd.DataFrame) for facet violin plots

    Parameters
    ----------
    tidy_df

    Returns
    -------

    """
	tidy_copy = copy.deepcopy(tidy_df)
	tidy_copy['stimname_row'] = tidy_df['stimname_col']
	tidy_copy['stimname_col'] = tidy_df['stimname_row']
	tidy_copy['stim_idx_row'] = tidy_df['stim_idx_col']
	tidy_copy['stim_idx_col'] = tidy_df['stim_idx_row']
	tidy_copy = tidy_copy.loc[tidy_copy['stimname_row'] != tidy_copy['stimname_col'], :]

	return tidy_df.append(tidy_copy)
#%%

def fix_rear_type(s):
	if 'pfo' in s:
		return 'pfo'
	else:
		return '1-6ol+EP'


# statfile.npz contains
"""
Load current datasets and filepaths from stat_files.npz

It includes:  

	pfo_datasets : namedtuple 
	odor_datasets : namedtuple
	stat_file_list : nesting listing of filepaths to suite2p/**/stat.npy
	flat_file_list : flat/unnested version of stat_file_list
"""
statfile_npz = np.load('stat_files.npz', allow_pickle=True)
print(statfile_npz.files)
flat_file_list = statfile_npz['flat_file_list']

tidy_corr = pd.DataFrame()
for i, f in enumerate(flat_file_list):
	df = feather.read_feather(f.with_name('tidy_corr_all_tril.feather'))
	# df['file_idx'] = i
	tidy_corr = tidy_corr.append(df)

tidy_corr['chunk'] = tidy_corr['chunk'].fillna(-1)

# make rear_type = 'pfo' or '1-6ol+EP'
rear_type = tidy_corr['rear_type'].apply(fix_rear_type)
tidy_corr['rear_type'] = rear_type
#%% ALL CONCENTRATIONS, POOLED
df_plot = tidy_corr.loc[(tidy_corr.corrtype=='cosine') & (tidy_corr.cellpop=='all'), :]
df_plot = reflect_tidy_corr(df_plot)

# Plot w/ all data, pooled across flies (all concentrations dumped in together)
g = sns.catplot(data=df_plot, x="stimname_col", y="value", hue="rear_type",
			                col="stimname_row", col_wrap=3,
			                palette="deep", legend_out=True, dodge=True,
                            height=4, aspect=0.75)
g.fig.suptitle("Cosine distance between odor pairs" +
               "\n8 flies, 26 movies, ~15-18 stimuli per movie, 1 point per trial-trial comparison"
               "\nconcentrations pooled (1e-05, 1e-04, 1e-03)")

plt.subplots_adjust(top=0.9, left=0.05, bottom=0.1)  # adjust spacing
			g.set_axis_labels("", 'cosine dist.')
			g.set(ylim=(-0.2, 1.0))
			g.map(plt.axhline, y=0, linewidth=0.5, color='k')
			#pdf.savefig(g.fig)
			#plt.show()
plt.show()
g.fig.savefig(Path.cwd().joinpath('doc', 'pooled_cosine_dist_pooled.png'))
g.fig.savefig(Path.cwd().joinpath('doc', 'pooled_cosine_dist_pooled.pdf'))
#%%
df_plot = tidy_corr.loc[(tidy_corr.corrtype=='cosine') & (tidy_corr.cellpop=='all') & (tidy_corr.conc==-5), :]
df_plot = reflect_tidy_corr(df_plot)


g = sns.catplot(data=df_plot, x="stimname_col", y="value", hue="rear_type",
			                col="stimname_row", col_wrap=3,
			                palette="pastel", legend_out=True, dodge=True,
                            height=4, aspect=0.75)
g.fig.suptitle("Cosine distance between odor pairs" +
               "\n odor concentrations = 1e-05" +
               "\n1 point per trial-trial comparison")

plt.subplots_adjust(top=0.9, left=0.05, bottom=0.1)  # adjust spacing
			g.set_axis_labels("", 'cosine dist.')
			g.set(ylim=(-0.2, 1.0))
			g.map(plt.axhline, y=0, linewidth=0.5, color='k')
			#pdf.savefig(g.fig)
			#plt.show()
plt.show()
g.fig.savefig(Path.cwd().joinpath('doc', 'pooled_cosine_dist_1e-05.png'))
g.fig.savefig(Path.cwd().joinpath('doc', 'pooled_cosine_dist_1e-05.pdf'))
#%%
df_plot = tidy_corr.loc[(tidy_corr.corrtype=='cosine') & (tidy_corr.cellpop=='all') & (tidy_corr.conc==-4), :]
df_plot = reflect_tidy_corr(df_plot)

# Plot w/ all data, pooled across flies (all concentrations dumped in together)
g = sns.catplot(data=df_plot, x="stimname_col", y="value", hue="rear_type",
			                col="stimname_row", col_wrap=3,
			                legend_out=True, dodge=True,
                            height=4, aspect=0.75)
g.fig.suptitle("Cosine distance between odor pairs" +
               "\n odor concentrations = 1e-04" +
               "\n1 point per trial-trial comparison")

plt.subplots_adjust(top=0.9, left=0.05, bottom=0.1)  # adjust spacing
			g.set_axis_labels("", 'cosine dist.')
			g.set(ylim=(-0.2, 1.0))
			g.map(plt.axhline, y=0, linewidth=0.5, color='k')
			#pdf.savefig(g.fig)
			#plt.show()
plt.show()
g.fig.savefig(Path.cwd().joinpath('doc', 'pooled_cosine_dist_1e-04.png'))
g.fig.savefig(Path.cwd().joinpath('doc', 'pooled_cosine_dist_1e-04.pdf'))
#%%
#%%
df_plot = tidy_corr.loc[(tidy_corr.corrtype=='cosine') & (tidy_corr.cellpop=='all') & (tidy_corr.conc==-3), :]
df_plot = reflect_tidy_corr(df_plot)

# Plot w/ all data, pooled across flies (all concentrations dumped in together)
g = sns.catplot(data=df_plot, x="stimname_col", y="value", hue="rear_type",
			                col="stimname_row", col_wrap=3,
			                palette="deep", legend_out=True, dodge=True,
                            height=4, aspect=0.75)
g.fig.suptitle("Cosine distance between odor pairs" +
               "\n odor concentrations = 1e-03" +
               "\n1 point per trial-trial comparison")

plt.subplots_adjust(top=0.9, left=0.05, bottom=0.1)  # adjust spacing
			g.set_axis_labels("", 'cosine dist.')
			g.set(ylim=(-0.2, 1.0))
			g.map(plt.axhline, y=0, linewidth=0.5, color='k')
			#pdf.savefig(g.fig)
			#plt.show()
plt.show()
g.fig.savefig(Path.cwd().joinpath('doc', 'pooled_cosine_dist_1e-03.png'))
g.fig.savefig(Path.cwd().joinpath('doc', 'pooled_cosine_dist_1e-03.pdf'))

#%%
grouped = tidy_corr.groupby(by=['corrtype', 'cellpop'])
#%%
df = grouped.get_group(('cosine', 'all'))
df.columns.to_list().remove(['stim_idx_row', 'stim_idx_col']

trial_grouper = df.groupby(by=['stimname_row', 'stimname_col', 'value',
        'flydate', 'flynum', 'mov', 'rear_type', 'chunk', 'conc']).mean()

	df.set_index(['flydate', 'flynum', 'mov', 'chunk', 'conc', 'stimname_col','stimname_row', 'rear_type'],
	             inplace=True)
	df.groupby(by=['corrtype', 'cellpop']).mean()

	tidy_corr_triu = tidy_corr_tril.append()
	tidy_corr_tril = tidy_corr_tril.append(feather.read_feather(f.with_name('tidy_corr_all_tril.feather')))

tidy_corr = pd.concat([tidy_corr_tril, tidy_corr_triu], ignore_index=True, verify_integrity=True)
#%%


# odor pair column, for seaborn plotting
#tidy_corr['odor_pair'] = [f"{row.stimname_row},{row.stimname_col}" for row in tidy_corr.itertuples()]
# %%
grouped = tidy_corr.groupby(by=['flydate', 'flynum', 'mov', 'chunk', 'conc', 'rear_type',
                                          'stimname_row', 'stimname_col', 'corrtype'],
                                      dropna=True, as_index=False)
tidy_corr_meanmov = grouped.mean()
tidy_corr_meanmov.dropna(inplace=True)

#%%
df_plot = tidy_corr_meanmov.loc[(tidy_corr_meanmov['corrtype'] == 'cosine')]
df_summary = df_plot.groupby(by=['odor_pair', 'stimname_row', 'stimname_col',
                                 'rear_type', 'conc']).mean()
df_summary.reset_index()


# %%

# %%
# corrtype, rear_type, conc, cellpop
df_plot = tidy_corr_meanmov.loc[(tidy_corr_meanmov['corrtype'] == 'pearson')]

g = sns.catplot(data=df_plot.sort_values(by=['stimname_row', 'stimname_col']),
                x="stimname_col", y="value",  col='stimname_row', col_wrap=3,
                 hue="conc", palette='muted')
plt.tight_layout(h_pad=2)
plt.show()
#%%
