library(here)
library(tidyverse)
library(ggpubr)
library(multcompView)

# Load accuracy lists
load(file = here("output", "accGP_Reml.RData"))
load(file = here("output", "accGP_LK.RData"))
load(file = here("output", "accGP_GK.RData"))

# Stacking lists
accGlobal <- c(accGP_Reml, accGP_LK, accGP_GK)

# A 2-way ANOVA with model type and modeling approach as different factors
# should be interesting

#  Data frame with all the information consolidated
acc <- data.frame(
  model = as.factor(rep(rep(str_remove(names(accGlobal)[1:5], "acc"), each = 10), 3)),
  approach = as.factor(rep(c("Reml", "LK", "GK"), each = 50)),
  CV_rep = as.factor(rep(1:10, length(accGlobal))),
  accuracy = unlist(accGlobal, use.names = FALSE)
)

# In order to assess whether the models differ significantly from each other, 
# we can employ an analysis of variance (ANOVA), using the CV repetition
# as a grouping factor (specified as a stratum via the Error(.) term), and
# model and approach as factors that may interact

modAOV <- aov(accuracy ~ model*approach + CV_rep, acc)
summary(modAOV)

# The interaction term is significant, implying that the model + approach
# combinations differ from each other, hence we can assess their pairwise
# differences via Tukey HSD, which will be illustrated in the last plot

#-------------------------------------------------------------#

# Color-blind friendly colors:
# Okabe-Ito palette (the standard reference palette for CVD-safety)
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")

#-----------------------------------------------------------------#
# Plot illustrating the indirect selection models relative to the 
# baseline field model, grouped by approach:

# Number of CV repetitions per (model + approach) combination
nRep <- length(unique(acc$CV_rep))

# Filtering for only the baseline Field model across all 3 frameworks
# As it will be used for horizontal contrast lines for each approach
Field_models <- acc |>
  filter(model == "Field") |>
  group_by(approach) |>
  summarise(mean_acc = mean(accuracy), se = sd(accuracy)/sqrt(nRep), 
            .groups = "drop")

# Filtering for remaining models for bar plots
IS_models <- acc |>
  filter(model != "Field") |>
  group_by(model, approach) |>
  summarise(mean_acc = mean(accuracy), se = sd(accuracy)/sqrt(nRep), 
            .groups = "drop")

ggplot(IS_models, aes(x = factor(model, levels = c("ColeoIS", "MesoIS", "MT_IS", "Idx")), 
                      y = mean_acc, 
                      fill = factor(model, levels = c("ColeoIS", "MesoIS", "MT_IS", "Idx")))) +
  geom_rect(
    data = Field_models,
    aes(ymin = mean_acc - se, ymax = mean_acc + se, xmin = -Inf, xmax = Inf),
    inherit.aes = FALSE, fill = "grey70", alpha = 0.4
  ) +
  geom_hline(
    data = Field_models,
    aes(yintercept = mean_acc),
    inherit.aes = FALSE, linetype = "solid", color = "black"
  ) +
  geom_text(
    data = Field_models,
    aes(x = -Inf, y = mean_acc, label = "Field"),
    inherit.aes = FALSE, hjust = -0.1, vjust = -0.5,
    size = 3.2, color = "black", fontface = "italic"
  ) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_acc - se, ymax = mean_acc + se), width = 0.2) +
  facet_wrap(~approach) +
  labs(x = NULL, y = "Predictive ability", fill = "Model type") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "grey85", color = "grey30", linewidth = 0.6),
    strip.text = element_text(face = "bold")
    ) +
  scale_fill_manual(values = okabe_ito, name = "Model type") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

#--------------------------------------------------------------------#

# Heatmap (5 x 3 grid comparing PAs across models and approaches)
summary_acc <- acc |>
  group_by(model, approach) |>
  summarise(mean_acc = mean(accuracy), se = sd(accuracy)/sqrt(nRep), 
            .groups = "drop")

ggplot(summary_acc, aes(x = approach, 
                        y = factor(model, levels = c("ColeoIS", "MesoIS", "MT_IS", "Idx", "Field")), 
                        fill = mean_acc)) +
  geom_tile() +
  geom_text(aes(label = round(mean_acc, 3)), color = "black") +
  scale_fill_viridis_c() +
  labs(x = "Approach", y = "Model", fill = "Mean Predictive Ability")

#--------------------------------------------------------------------------#

# "Delta" plot comparing each approach's IS models to their respective
# baseline field model, computing the differences for each CV repetition
acc_delta <- acc |>
  group_by(approach, CV_rep) |>
  mutate(Field_acc = accuracy[model == "Field"]) |>
  ungroup() |>
  filter(model != "Field") |>
  mutate(delta_acc = accuracy - Field_acc)

# Summarizing for plotting
delta_summary <- acc_delta |>
  group_by(model, approach) |>
  summarise(mean_delta = mean(delta_acc), se = sd(delta_acc)/sqrt(nRep), 
            .groups = "drop")

ggplot(delta_summary, aes(x = factor(model, levels =  c("ColeoIS", "MesoIS", "MT_IS", "Idx")), 
                          y = mean_delta,
                          fill = factor(model, 
                          levels = c("ColeoIS", "MesoIS", "MT_IS", "Idx")))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_delta - se, ymax = mean_delta + se), width = 0.2) +
  facet_wrap(~approach) +
  scale_fill_manual(values = okabe_ito, name = "Model type") +
  labs(x = NULL, y = "PA loss") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "grey85", color = "grey30", linewidth = 0.6),
    strip.text = element_text(face = "bold")
  )

#------------------------------------------------------------------------#
# 15-model ranking - lollipop plot

# In order to conduct pairwise comparisons, we can conduct the Tukey's HSD 
# (Honestly Significant Difference):
(tukeyMod <- TukeyHSD(modAOV, which = "model:approach"))

# Adding CLD (Compact Letter Display) to the lollipop plot:
tukey_pvals <- tukeyMod$model[, "p adj"]
lettersComp <- multcompLetters(tukey_pvals)

# Single column in the original acc dataset combining model and approach
# with ":" as a separator to ensure compatibility with lettersComp (via y_pos)
acc <- acc |>
  mutate(mod_app = as.factor(paste(model, approach, sep = ":")))

# Data frame with letter positions:
letterDF <- data.frame(
  # Letters placed right above the max accuracy value for each model
  mod_app = names(lettersComp$Letters),
  letter = lettersComp$Letters
)
rownames(letterDF) = NULL

# Data frame with mean accuracy values for each model + approach
rank_summary <- acc |>
  # Redundant simply because we want to keep all columns
  group_by(model, approach, mod_app) |>
  summarise(mean_acc = mean(accuracy), se = sd(accuracy)/sqrt(nRep), 
            .groups = "drop") |>
  # Merging rank_summary to letterDF
  left_join(letterDF, by = "mod_app") |>
  arrange(mean_acc)
  

# Modifying the model + accuracy column for better display in the lollipop plot
rank_summary <- rank_summary |>
                mutate(mod_app = str_replace_all(mod_app, ":", " \u2013 "))

# To ensure that the lines in the plot will be in ascending order of 
# predictive ability
levels_mod_app <- levels(reorder(rank_summary$mod_app, rank_summary$mean_acc))
  
ggplot(rank_summary, aes(x = mean_acc, y = factor(mod_app, levels = levels_mod_app), 
                                                  color = approach)) +
  geom_segment(aes(x = 0.4, xend = mean_acc, y = mod_app, yend = mod_app), 
               linewidth = 0.6) +
  geom_point(size = 3.5) +
  geom_errorbarh(aes(xmin = mean_acc - se, xmax = mean_acc + se), height = 0.2) +
  geom_text(aes(x = mean_acc + se + 0.01, label = letter), hjust = 0, size = 3.2,
            show.legend = FALSE) +
  scale_color_manual(values = c("#1b9e77", "#d95f02","#7570b3"),
                     name = "Approach") +
  labs(x = "Predictive ability", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

