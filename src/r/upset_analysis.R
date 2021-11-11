# ---
# title: "upset responders"
# author: remy yang
# date: august 10, 2021
# output: html_document
# ---


# imports -----------------------------------------------------------------
#+ results=FALSE
library(tidyverse)
library(ggplot2)
library(ComplexUpset)
library(ggplot2movies)
library(reticulate)
library(yaml)
library(R.utils)
library(ggsci)

# load .npy files ---------------------------------------------------------
reticulate::use_condaenv('suite2p')
np <- import("numpy")

# ComplexUpset plot: Counts -------------------------------------------------------
upset_counts <- function(dat){
    K <- nrow(dat)
    ComplexUpset::upset(
        dat, colnames(dat), name='responders', width_ratio=0.3, min_size=0, keep_empty_groups=TRUE, wrap=TRUE,
        sort_sets=FALSE,
        sort_intersections_by=c('degree'), sort_intersections='ascending',
        group_by='sets',
        matrix=intersection_matrix(geom=geom_point(size=2)),
        queries=list(upset_query(group='EP', color='blue', fill='blue'),
                     upset_query(group='1-6ol', color='orange', fill="orange"),
                     upset_query(group='1-6ol+EP',color='green', fill="green")),
        base_annotations=list(
            #'exclusive\nintersection'=intersection_size(mode='distinct', text=list(size=2.5, angle=90, vjust=0.5, hjust=-0.15)),
            'inclusive\nintersection'=intersection_size(mode='intersect', text=list(size=2.5, angle=90, vjust=0.5, hjust=-0.15)),
            'inclusive\nunion'=intersection_size(mode='union', text=list(size=2.5, angle=90, vjust=0.5, hjust=-0.15)),
            'Intersection\nratio'=intersection_ratio(mode='intersect',
                                                     denominator_mode = "union",
                                                     text_mapping=aes(label=!!upset_text_percentage(mode='intersect', digits=1)),
                                                     text=list(size=2.5, angle=90, vjust=0.5, hjust=-0.15)
            )
        ),
        set_sizes=(
            upset_set_size() +
            geom_text(aes(label=scales::percent(round((..count..)/sum(..count..),3))), hjust=1.1, stat='count', size=2.5) +
            expand_limits(y=K) +
            theme(axis.ticks.x=element_line())
        )
    )
}

upset_percentages <- function(dat){
    K <- nrow(dat)
    rng <- NULL
    ComplexUpset::upset(
        dat, colnames(dat), name='responders', width_ratio=0.3, min_size=0, keep_empty_groups=TRUE, wrap=TRUE,
        sort_sets=FALSE,
        group_by='sets',
        matrix=intersection_matrix(geom=geom_point(size=2)),
        queries=list(upset_query(group='EP', color='orange', fill='orange'),
                     upset_query(group='1-6ol', color='blue', fill="blue"),
                     upset_query(group='1-6ol+EP',color='green',  fill="green")),
        sort_intersections_by=c('degree'),  sort_intersections='ascending',
        base_annotations=list(
            'exclusive\nintersection'=intersection_size(mode='distinct',
                                                        text=list(size=2.5, angle=90, vjust=0.5, hjust=-0.15),
                                                        text_mapping=aes(label=paste0(round(!!get_size_mode('exclusive_intersection') * (100/K),1), '%')))+
            scale_y_continuous(limits=rng,breaks=c(0,.25,.5,.75,1)*K, labels=scales::label_percent(scale=(100/K))),
            'inclusive\nintersection'=intersection_size(mode='intersect',
                                                        text=list(size=2.5, angle=90, vjust=0.5, hjust=-0.15),
                                                        text_mapping=aes(label=paste0(round(!!get_size_mode('intersect') * (100/K),1), '%')))+
                scale_y_continuous(limits=rng,breaks=c(0,.25,.5,.75,1)*K, labels=scales::label_percent(scale=(100/K))),
            'inclusive\nunion'=intersection_size(mode='union',
                                                 text=list(size=2.5, angle=90, vjust=0.5, hjust=-0.15),
                                                 text_mapping=aes(label=paste0(round(!!get_size_mode('inclusive_union') * (100/K),1), '%'))) +
            scale_y_continuous(limits=rng,breaks=c(0,.25,.5,.75,1)*K, labels=scales::label_percent(scale=(100/K))),
            'Intersection\nratio'=intersection_ratio(mode='intersect',
                                                     denominator_mode = 'union',
                                                     text=list(size=2.5, angle=90, vjust=0.5, hjust=-0.1),
                                                     text_mapping=aes(label=!!upset_text_percentage(mode='intersect', digits=1))
            )
        ),
        set_sizes=(
            upset_set_size()
            + geom_text(aes(label=scales::percent(round((..count..)/K,3))), hjust=1.1, stat='count', size=3)
            + expand_limits(y=100)
            + theme(axis.ticks.x=element_line())
        ),
    )
}

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

    # upset plot data table
    #df_tuning <- py_to_r(binresp$df_tuning_logical) %>%

    df_tuning <- as.data.frame((py_to_r(binresp$df_tuning) > 1))*1




    # subtitle string
    s = sprintf("date:%s | fly #:%s | movie:%s | reared:%s | concentration:%d", fly_date, fly_num, movie, reared, params[["movie_conc"]])

    cap = "params:"
    for (item in names(params)){
        cap <- paste(cap, "\n",  item, ": ",toString(jsonlite::toJSON(params[[item]],auto_unbox=TRUE)))
    }
    # make plots --------------------------------------------------------------

    # g1: count plots + pfo
    g1 <- upset_counts(df_tuning) +
        labs(title="Responders - intersections by type (counts)",
             subtitle=s,
             caption=cap) +
        theme(plot.caption = element_text(hjust=0))

    # g2:  % plot + pfo
    g2 <- upset_percentages(df_tuning) +
        labs(title="Responders - intersections by type (%)",
             subtitle=s,
             caption=cap) +
        theme(plot.caption = element_text(hjust=0))

    # g3: count plots, no pfo
    g3 <- upset_counts(select(df_tuning, !"pfo")) +
        labs(title="Responders - intersections by type (counts), no pfo",
             subtitle=s,
             caption=cap) +
        theme(plot.caption = element_text(hjust=0))


    # g4: make % plot, no pfo
    g4 <-  upset_percentages(select(df_tuning, !"pfo")) +
        labs(title = "Responders - intersections by type (%), no pfo",
             subtitle = s,
             caption=cap) +
        theme(plot.caption = element_text(hjust=0))

    # save plots --------------------------------------------------------------

    save_individual_plots=FALSE
    save_multipage_pdf=TRUE

    if (save_individual_plots){
        w_in <-  11
        h_in <-  8
        ggsave(file.path(savesubdir, 'upset_counts.png'), plot=g1, width=w_in, height=h_in, units="in")
        ggsave(file.path(savesubdir, 'upset_counts.pdf'), plot=g1, width=w_in, height=h_in, units="in")

        ggsave(file.path(savesubdir, 'upset_percentages.png'), plot=g2, width=w_in, height=h_in, units="in")
        ggsave(file.path(savesubdir, 'upset_percentages.eps'), plot=g2, width=w_in, height=h_in, units="in")

        ggsave(file.path(savesubdir, 'upset_counts_nopfo.png'), plot=g3, width=w_in, height=h_in, units="in", device="png")
        ggsave(file.path(savesubdir, 'upset_counts_nopfo.pdf'), plot=g3, width=w_in, height=h_in, units="in", device="pdf")

        ggsave(file.path(savesubdir, 'upset_percentages_nopfo.png'), plot=g4, width=w_in, height=h_in, units="in", device="png")
        ggsave(file.path(savesubdir, 'upset_percentages_nopfo.pdf'), plot=g4, width=w_in, height=h_in, units="in", device="pdf")
    }

    if (save_multipage_pdf){
        pdf(file.path(savesubdir, 'upset_responders_groupedsets.pdf'),width=11, height=8)
        print(g1)
        print(g2)
        print(g3)
        print(g4)
        dev.off()

    }
}



# make count plot + pfo ---------------------------------------------------

#text = list(size=3)
#)  + scale_y_continuous(labels=scales::percent_format(scale=1)),


# R example - movies
# ComplexUpset::upset(
#     movies, genres, name='genre', width_ratio=0.1, min_size=10,
#     mode='inclusive_union',
#     base_annotations=list(
#         # with manual aes specification:
#         'Intersection size'=intersection_size(text_mapping=aes(label=paste0(round(
#                                             !!get_size_mode('inclusive_intersection') * (100/K)
#                                             ), '%'))),
#         # using shorthand:
#         'Intersection ratio'=intersection_ratio(text_mapping=aes(label=!!upset_text_percentage()))
#     )
