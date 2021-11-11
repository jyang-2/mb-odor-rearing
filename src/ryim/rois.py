"""
This module contains functions for drawing rois, and writing tiff label images, using suite2p outputs
"""


#%%
folder = r"G:\HongLab @ Caltech Dropbox\Remy\mb_odor_rearing\data\processed_data\2021-08-30\1\movie_002\1\suite2p"
folder = Path(folder)
print(f"suite2p folder: {folder}")
from natsort import natsorted, realsorted


statfile_list = realsorted(list(folder.rglob("**/plane*/stat.npy")))
for f in statfile_list:
    print(f.relative_to(Path.cwd()))
#%%
statfile_list = [str(item) for item in statfile_list]
#%%
statfile_list = list(map(lambda x: Path(x), realsorted(statfile_list)))
#%%
stat_list = []
for f in statfile_list:
    stat_list.append(np.load(f, allow_pickle=True))
#%%
stat_list = list(map(lambda x: np.load(x, allow_pickle=True), statfile_list))
ops_list = list(map(lambda x: np.load(x.with_name('ops.npy'), allow_pickle=True).item(), statfile_list))
#%%
roi_3D_array = stats.ROI.stats_dicts_to_3d_array(stat_list[0], 256, 256, label_id=True)
#%%
import matplotlib.pyplot as plt




def draw_single_roi(ops, stat):
    im = np.zeros((ops['Ly'], ops['Lx']))
    ypix = stat['ypix'][~stat['overlap']]
    xpix = stat['xpix'][~stat['overlap']]
    im[ypix, xpix] = stat['lam']
    return im

label_img = draw_single_roi(ops_list[plane_idx], stat_list[plane_idx][0])
plt.imshow(label_img)
plt.show()

#%%


def roilist_to_3d_array(ops, stat):
    stack = np.zeros((len(stat), ops['Ly'], ops['Lx']))
    for z, stat0 in enumerate(stat):
        stack[z, :, :] = draw_single_roi(ops, stat0)
    return stack


label_stack = roilist_to_3d_array(ops_list[plane_idx], stat_list[plane_idx])
#%%
def draw_lambda_roi(ops, stat, ncells):
    im = np.zeros((ops['Ly'], ops['Lx']))
    for n in range(0, ncells):
        ypix = stat[n]['ypix'][~stat[n]['overlap']]
        xpix = stat[n]['xpix'][~stat[n]['overlap']]
        im[ypix, xpix] = stat[n]['lam']
    return im

def draw_label_roi(ops, stat, ncells):
    im = np.zeros((ops['Ly'], ops['Lx']), dtype=np.int)
    for n in range(0, ncells):
        ypix = stat[n]['ypix'][~stat[n]['overlap']]
        xpix = stat[n]['xpix'][~stat[n]['overlap']]
        im[ypix, xpix] = n + 1
    return im

def draw_plane_rois(ops, stat, roitype="lambda"):
    ncells = len(stat)
    if roitype == 'lambda':
        im = draw_lambda_roi(ops, stat, ncells)
    else:
        im = draw_label_roi(ops, stat, ncells)
    return im

plane_idx = 0
label_img = draw_plane_rois(ops_list[plane_idx], stat_list[plane_idx], roitype="label")
lambda_img = draw_plane_rois(ops_list[plane_idx], stat_list[plane_idx], roitype="lambda")
#%%
from skimage import exposure

fig, axs = plt.subplots(1, 2, figsize=(16,9))
alpha = colors.Normalize(0, 0.8, clip=True)(lambda_img)
axs[0].pcolormesh(label_img, cmap=cmap, alpha=alpha)
axs[1].pcolormesh(lambda_img, cmap='gray')
plt.show()

#%% plot roi overlay
import matplotlib.colors as colors



n_cells = len(stat_list[plane_idx])
#%% make colormap for label -> rgba mapping
cmap = plt.cm.get_cmap('viridis', n_cells + 1)
cmap = cmap(np.linspace(np.arange(n_cells+1)))
cmap = colors.ListedColormap()
cmap = cmap.colors
rng = np.random.default_rng()
rng.shuffle(cmap)
cmap[0,:]=[1, 1, 1, 0]
cmap = colors.ListedColormap(cmap)
cmap.set_under(color='w', alpha=0)
#%%
color1 = label2rgb(label_img, image=None, bg_label=0, bg_color=None, kind='avg')


meanImg = ops_list[plane_idx]['meanImg']

alphas = np.clip(alphas, 0.0, 0.5)  # alpha value clipped at the bottom at .4
#%%
#alpha_img =colors.Normalize(0, 0.8, clip=True)(lambda_img)
alpha_img = np.clip(lambda_img/2, 0.0, 0.3)  # alpha value clipped at the bottom at .4
#%%
colorImg = np.take(cmap.colors, label_img, axis=0)
colorImg[:,:,3] = lambda_img/3
#%%
from skimage import exposure
# Create the figure and image
# Note that the absolute values may be slightly different
img_adapteq = exposure.equalize_adapthist(colors.Normalize()(meanImg), clip_limit=0.03)
#%%
fig, ax = plt.subplots()
ax.imshow(img_adapteq, cmap='gray')
ax.imshow(colorImg)
#ax.imshow(color1, alpha=alpha)
plt.show()
#%% shuffled colormap
vals = np.linspace(0,1,256)
np.random.shuffle(vals)
