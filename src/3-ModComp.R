library(here)
library(asreml)
library(dplyr)
library(tidyr)
library(ggplot2)
library(asremlPlus)


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

lapply(accs_List, function(a) {
  qqnorm(a)
  qqline(a)
})
# Normality seems like a reasonable assumption


# Model mean accuracies and standard deviations for these
accs <- data.frame(
  Model = names(accs_List),
  Accuracy = unlist(lapply(accs_List, mean)),
  StdDev = unlist(lapply(accs_List, sd))
)

# Bar chart
ggplot(accs, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_col() +
  theme_minimal() +
  labs(title = "Mean Accuracies", y = "Accuracy", x = "") +
  geom_errorbar(aes(ymin = Accuracy - StdDev, ymax = Accuracy + StdDev)) +
  theme(axis.text.x = element_blank()) +
  geom_text(aes(label = round(Accuracy, 2), vjust = -0.7)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# In order to compare the model accuracies, we will perform analysis of variance,
# and Tukey HSD for pairwise differences

accDF <- as.data.frame(accs_List)

# Adding a column to represent the 5-fold CV repetition that groups 
# the accuracy values across models
accDF$repID <- row.names(accDF)

# Converting data frame to long format to prepare for aov
accDF <- accDF |>
  pivot_longer(-repID, names_to = "model", values_to = "accuracy") |>
  mutate(model = as.factor(model), repID = as.factor(repID))

# Accuracy is modeled in terms of the model plus a fluctuation due to the
# specific 5-fold CV repetition
modAOV <- aov(accuracy ~ model + Error(repID), accDF)
summary(modAOV)

# Since "model" has a significant effect, there is evidence that at least 
# two models differ significantly

# Note: HSD corresponds to "honestly significant difference"
# For TukeyHSD, we fit a mixed model to the accuracy ~ model relationship, so
# repID can be treated as a random effect (the Error notation from aov has poor
# compatibility with TukeyHSD)


modLMM <- asreml(
  fixed = accuracy ~ model,
  random = ~ factor(repID), 
  data = accDF
)

pred_obj <- predictPlus(
  # Differences across models
  classify = "model",
  asreml.obj = modLMM,
  # Wald tests whether the differences are statistically equal to 0
  wald.tab = wald.asreml(modLMM, denDF = "numeric")$Wald
)

# Note: LSD stands for "least significant differences"
# Pairwise comparisons
pred_obj$differences    # LSD-based differences
pred_obj$LSD            # LSD threshold

# The assigned LSD is 0.005948599 (as seen in pred_obj$LSD), and all values
# in pred_obj$differences are larger than this value (in module), hence
# all the models differ significantly from each other

sum(abs(pred_obj$differences) < 0.005948599)
# There are 6 differences below the threshold, but they correspond to
# the models compared to themselves

# Now, to rank the models, we can add up how many negative differences are
# found in each row, which correspond to how many models perform better than
# the one represented in that given row (that + 1 is their rank)

rowSums(pred_obj$differences < 0) + 1 # model ranking
