# Title     : DoOR 2.0.1 - odor map analysis
# Objective : Identify odor set for rearing experiment
# Created by: remy
# Created on: 7/19/21

library(ggplot2)
library(ggrepel)
library(patchwork)
library(DoOR.data)
library(DoOR.functions)
library(embed)
library(factoextra)


load_door_data(nointeraction = TRUE)

hallem_df <- get_dataset("Hallem.2006.EN", na.rm=TRUE)
#                 door_response_matrix - globally normalized
# door_response_matrix_non_normalized  - is a version of the consensus data that is not globally normalized,
#                                        meaning that responses are scaled [0,1] within each responding unit
#                                        (receptor, sensory neuron, glomerulus…).

odors <- trans_id(c("pentyl acetate", "ethyl propionate", "1-hexanol", "1-pentanol", "pentanoic acid"),
           from = "Name")
dplot_across_ru("ethyl propionate", tag = "Name", base_size = 10)
dplot_across_ru("ethyl propionate", tag = "Name", base_size = 10)


# AL map
p1 <- dplot_al_map(trans_id("1-hexanol", from="Name"), base_size = 10)
p2 <- dplot_al_map(trans_id("1-pentanol", from="Name"), base_size = 10)
p3 <- dplot_al_map(trans_id("ethyl propionate", from="Name"), base_size = 10)
p <- p1 / p2 /p3


# dplot_tuningCurve() can as well be used to visualize the ensemble of responding units that is activated by a given odorant. Therefore, we specify an odorant name instead of a response unit name:
t1 <- dplot_tuningCurve(odorant = trans_id("1-hexanol", from="Name"), base_size = 10)
t2 <- dplot_tuningCurve(odorant = trans_id("ethyl propionate", from="Name"), base_size = 10)
t <- t1 / t2

#
library(tidyverse)
library(tidymodels)

pca_rec <- recipe(~., data = hallem_df) %>%
          update_role(colnames(hallem_data)[1:5], new_role = "id") %>%
          step_normalize(all_predictors()) %>%
          step_pca(all_predictors())
pca_prep <- prep(pca_rec)



# plot pca component vectors
tidied_pca <- tidy(pca_prep, 2)
tidied_pca %>%
  filter(component %in% paste0("PC", 1:5)) %>%
  mutate(component = fct_inorder(component)) %>%
  ggplot(aes(value, terms, fill = terms)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~component, nrow = 1) +
  labs(y = NULL)

options(ggrepel.max.overlaps = Inf)

tidied_pca %>%
  filter(component %in% paste0("PC", 1:4)) %>%
  group_by(component) %>%
  top_n(8, abs(value)) %>%
  ungroup() %>%
  mutate(terms = reorder_within(terms, abs(value), component)) %>%
  ggplot(aes(abs(value), terms, fill = value > 0)) +
  geom_col() +
  facet_wrap(~component, scales = "free_y") +
  scale_y_reordered() +
  labs(
    x = "Absolute value of contribution",
    y = NULL, fill = "Positive?"
  )


pca_df = juice(pca_prep)
in_odors = c("1-hexanol", "1-pentanol", "ethyl propionate", "benzaldehyde")

# family = "IBMPlexSans"
pca_plot <- juice(pca_prep) %>%
  ggplot(aes(PC1, PC2, label = Name)) +
  #geom_text_repel() +
  geom_point(aes(color = Class), alpha = 0.7, size = 2) +
  #geom_label(data=hallem_df[aes(fill = factor(Class)), colour = "white", fontface = "bold") +
  geom_text(check_overlap = TRUE, hjust = "inward", ) +
  labs(color = NULL)

# ===============================
library(embed)

umap_rec <- recipe(~., data = cocktails_df) %>%
  update_role(name, category, new_role = "id") %>%
  step_normalize(all_predictors()) %>%
  step_umap(all_predictors())

umap_prep <- prep(umap_rec)

umap_prep
# -------------------------------
# run U-MAP
# -------------------------------
swissUmap <- select(hallem_df, -Status) %>%
             as.matrix() %>%
             umap(n_neighbors = 7, min_dist = 0.1,
                  metric = "manhattan", n_epochs = 200, verbose = TRUE)
