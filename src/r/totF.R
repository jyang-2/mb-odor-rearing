

# imports -----------------------------------------------------------------
#+ results=FALSE
library(tidyverse)
library(ggplot2)
library(ComplexUpset)
library(ggplot2movies)
library(reticulate)
library(yaml)
library(R.utils)
library(patchwork)

# load .npy files ---------------------------------------------------------
reticulate::use_condaenv('suite2p')
np <- import("numpy")

# run upsetplots ----------------------------------------------------------

# get list of respvecs.npy files
files = list.files(path="/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing",
                   pattern='respvecs.npy', recursive=TRUE, full.names=TRUE)

for(filepath in files){
  savesubdir = normalizePath(dirname(filepath))
  PROC_DIR = dirname(savesubdir)
  FLY_DIR = dirname(PROC_DIR)
  movie = basename(PROC_DIR)
  fly_num = basename(FLY_DIR)
  fly_date = basename(dirname(FLY_DIR))
  print(savesubdir)


  metadata <- yaml::read_yaml(file.path(FLY_DIR, 'metadata.yaml'))
  reared <- metadata[["Fly info"]][[fly_num]][["reared"]]
  df_stimulus <- read_csv(file.path(savesubdir, 'df_stimulus.csv'))
  params <- yaml::read_yaml(file.path(savesubdir, 'params.yaml'))
  responses <- np$load(file.path(savesubdir, 'respvecs.npy'), allow_pickle=TRUE)[[1]]
  binresp <- np$load(file.path(savesubdir, 'binresp.npy'), allow_pickle=TRUE)[[1]]
  traces <- np$load(file.path(savesubdir, 'traces.npy'), allow_pickle=TRUE)[[1]]
  trials_and_timing <- np$load(file.path(savesubdir, 'trials_and_timing.npy'), allow_pickle=TRUE)[[1]]

  df_frametimeinfo <- py_to_r(trials_and_timing$df_frametimeinfo)

  # total fluorescence
  totF <- colSums(traces$df)
  stimname <- ordered(df_stimulus$stimname, levels=unique(df_stimulus$stimname))
  s = sprintf("date:%s | fly #:%s | movie:%s | reared:%s | concentration:%d", fly_date, fly_num, movie, reared, params[["movie_conc"]])

  df_join <-  df_frametimeinfo %>%
    mutate(totF = totF) %>%
    mutate(rep = as.integer(trialtrace_idx %% 6)) %>%
    inner_join(df_stimulus, by=c("trialtrace_idx"="stim_idx")) %>%
    mutate(frame_times = round(frame_times, digits=1)) %>%
    mutate(stim_ict = round(stim_ict, digits=1)) %>%
    mutate(frame_times_zeroed = frame_times - stim_ict)


  # l1 <- ggplot(df_join, aes(x=frame_times_zeroed, y=totF, group=trialtrace_idx, color=stimname)) +
  #   geom_line() +
  #   facet_wrap(~ stimname, ncol = 5)
  # l1 <- l1 + labs(title=s, subtitle="totF (sum of df for all suite2p cells)", )
  #print(l1)

   l1 <- ggplot(df_join, aes(x=frame_times_zeroed, y=totF, color=stimname)) +
     geom_line(aes(group=trialtrace_idx), alpha=0.4)+
     geom_smooth(method="auto",span=0.3, se=TRUE, fullrange=TRUE, level=0.95) +
     facet_wrap(~ stimname, ncol=5)
   l1 <-l1 + labs(title=s, subtitle="totF (sum of df for all suite2p cells)", )
   #print(l2)

   l2 <- ggplot(df_join, aes(x=frame_times_zeroed, y=totF, color=stimname)) +
     #geom_point(alpha=0.1) +
     geom_line(aes(group=trialtrace_idx), alpha=0.5+
     geom_smooth(method="auto",span=0.3, se=TRUE, fullrange=TRUE, level=0.95)
   l2 <- l2 + labs(title=lapply(strwrap(s, width = 30, simplify = FALSE), paste, collapse="\n"),
          subtitle="totF (sum of df for all suite2p cells)", )


  save_individual_plots <- TRUE
  if (save_individual_plots){
    w_in <-  10
    h_in <-  3

    ggsave(file.path(savesubdir, 'tot_F.png'), plot=l1, width=w_in, height=h_in, units="in", device="png")
    ggsave(file.path(savesubdir, 'tot_F.pdf'), plot=l1, width=w_in, height=h_in, units="in", device="pdf")

    ggsave(file.path(savesubdir, 'tot_F_overlay.pdf'), plot=l2, width=5, height=5, units="in", device="pdf")
    ggsave(file.path(savesubdir, 'tot_F_overlay.png'), plot=l2, width=5, height=5, units="in", device="png")

  }
}