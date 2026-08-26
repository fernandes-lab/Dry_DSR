library(here)
library(tidyverse)
library(ggpubr)
library(multcompView)

###############################################################
##                Assessing model accuracies                 ##
###############################################################

# save(accs_List, file = here("output", "accs_List.RData"))

# The assessment uses output from this source file, from 2-GBLUP_Idx,
# and from 2-RKHS.R

load(file = here("output", "accs_List.RData"))

# Testing the accuracy vectors for normality to check whether t-tests are
# viable
lapply(accs_List, shapiro.test)

par(mfrow=c(2, 3))
lapply(accs_List, function(a) {
  qqnorm(a)
  qqline(a)
})
# Normality seems like a reasonable assumption

# Note: accuracy (predictive ability) in our context translates to
# ability to rank genotypes

#  Data frame with the respective models, accuracies and CV repetitions
acc_long <- data.frame(
  model = as.factor(rep(str_remove(names(accs_List), "acc"), each = 10)),
  accuracy = unlist(accs_List, use.names = FALSE),
  CV_rep = as.factor(rep(1:10, length(accs_List)))
)

# In order to assess whether the models differ significantly from each other, 
# we can employ an analysis of variance (ANOVA), using the CV repetition
# as a grouping factor (specified as a stratum via the Error(.) term):

modAOV <- aov(accuracy ~ model + CV_rep, acc_long)
summary(modAOV)

# We can see that model is a significant term to the model, implying that
# the models statistically differ from each other

# In order to conduct pairwise comparisons, we can conduct the Tukey's HSD 
# (Honestly Significant Difference):
(tukeyMod <- TukeyHSD(modAOV, which = "model"))

# Generate letter groupings for boxplot:
tukey_pvals <- tukeyMod$model[, "p adj"]
lettersComp <- multcompLetters(tukey_pvals)

# Data frame with letter positions:
y_pos = tapply(acc_long$accuracy, acc_long$model, max) + 0.02
letterDF <- data.frame(
  # Letters placed right above the maximum value for each
  # model 
  # max accuracy value across each model
  model = names(y_pos),
  letter = lettersComp$Letters,
  y.pos = y_pos
)
rownames(letterDF) = NULL
rm(y_pos)

# Boxplot for model predictive ability comparison:

ggboxplot(acc_long,
          x     = "model",
          y     = "accuracy",
          fill  = "model") +
  geom_text(data  = letterDF,
            aes(x = model, y = y.pos, label = letter),
            size  = 5,
            fontface = "bold") +
  labs(y    = "Predictive Ability",
       x    = "Model") +
  theme_classic() +
  theme(legend.position = "none")



