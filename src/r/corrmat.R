# Title     : TODO
# Objective : TODO
# Created by: remy
# Created on: 7/17/21

library(reticulate)
library(tidyverse)
library(reticulate)
library(readxl)
library(corrr)
library(RColorBrewer)
library(superheat)
library(ggplot2)
library(patchwork)

use_condaenv('suite2p', required = TRUE)

fname <- file.path(getwd(), "data", "processed_data", "2021-06-26", 1, "movie_001", "metadata.xlsx")
df_stimulus <- read_excel(fname, sheet="stimulus")


fname <- "/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/data/processed_data/2021-06-26/1/movie_001/__win30__prctile20/distances.npy"
np <- import("numpy")
dist <- np$load(fname, allow_pickle=TRUE)[[1]]
corr_cosine <- dist$cosine

responses <- np$load("/media/remy/remy-storage/Remy's Dropbox Folder/HongLab @ Caltech Dropbox/Remy/mb_odor_rearing/data/processed_data/2021-06-26/1/movie_001/__win30__prctile20/responses.npy",
                     allow_pickle=TRUE)[[1]]

corrmat = corrr::correlate(peakresp_amp=responses$peakresp_amp)

df_peakresp <- py_to_r(responses$df_peakresp)
df_corrmat <-

df <- correlate(df_peakresp) %>%
      shave() %>%
      stretch(na.rm=TRUE) %>%

df <- separate(x, c("stim_idx_a", "stimname_a"), sep="\\s") %>%
      separate(y, c("stim_idx_b", "stimname_b"), convert=TRUE, sep="\\s")

pattern <- "(\\d+)\\1\\s+([:alnum:])\\2"
df$x %>% str_subset(pattern) %>%
  str_match(pattern)

# # pal(100)
# superheatmap1 <- superheat(xs,
#           scale=FALSE,
#           row.dendrogram = TRUE,
#           heat.pal=pal(100),
#           grid.hline.col = "white",
#           grid.vline.col = "white",
#           bottom.label.text.angle = 90,
#           bottom.label.text.size = 4,
#           left.label.text.size = 4,
#           print.plot=TRUE)

clrs <-rev( RColorBrewer::brewer.pal(3, "RdBu"))
pal <- colorRampPalette(clrs)

peakresp_amp <- responses[["peakresp_amp"]]
stimname <- ordered(df_stimulus$stimname, levels=unique(df_stimulus$stimname))
#ord = order(stimname)
test[ , order(names(test))]

sc <- str_c(df_stimulus$stimname, ', ', df_stimulus$stim_idx)

df_peakresp_amp <- data.frame(peakresp_amp)
colnames(df_peakresp_amp) <- df_stimulus$stimname
df_peakresp_amp

rmat_pearson <- cor(df_peakresp_amp[, order(names(df_peakresp_amp))])

s1 <- superheat(as.matrix(correlate(peakresp_amp, use = "pairwise.complete.obs")),
          #order.rows=ord,
          #order.cols=ord,
          #scale=TRUE,
          #heat.lim = c(-0.8, 0.8),
          #pretty.order.cols=TRUE,
          #pretty.order.rows=TRUE,
          #row.dendrogram = TRUE,
          #col.dendrogram = TRUE,
          heat.pal=pal(100),
          grid.hline.col = "white",
          grid.vline.col = "white",
          bottom.label.text.angle = 90,
          bottom.label.text.size = 2,
          left.label.text.size = 2,
          membership.rows = stimname,
          membership.cols = stimname,
          print.plot=TRUE)

s2 <- superheat(rmat_pearson,
                #order.rows=ord,
                #order.cols=ord,
                scale=FALSE,
                heat.lim = c(-1, 1),
                pretty.order.cols=TRUE,
                pretty.order.rows=TRUE,
                #row.dendrogram = TRUE,
                #col.dendrogram = TRUE,
                heat.pal=pal(100),
                grid.hline.col = "white",
                grid.vline.col = "white",
                bottom.label.text.angle = 90,
                bottom.label.text.size = 2,
                left.label.text.size = 2,
                #membership.rows = stimname,
                #membership.cols = stimname,
                print.plot=TRUE)



corrgram(rmat_pearson, colors = c("red", "blue"), order=NULL, lower.panel=panel.shade,
         upper.panel=NULL, text.panel=panel.txt,
         main="Car Milage Data (unsorted)")

rplot(rmat_pearson, colors=c("red, blue"))

stimname <- ordered(df_stimulus$stimname, levels=unique(df_stimulus$stimname))

rplot(
  rmat_pearson,
  legend = TRUE,
  shape = 16,
  colours =clrs, #c("indianred2", "white", "skyblue1"),
  print_cor = TRUE,
)