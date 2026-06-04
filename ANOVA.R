# 2-way ANOVA testing of phenotypic traits

# 1. Check assumptions
shapiro.test(df_matsuda_index_female_full$matsuda_index) # Checking the normality assumption of a dataset
leveneTest(matsuda_index ~ Genotype * Diet, data = df_matsuda_index_female_full) # Checking the homogeneity of variance

# 2. Two-way ANOVA
model <- aov(matsuda_index ~ Genotype * Diet, data = df_matsuda_index_female_full)
summary(model)

# 3. If the interaction is significant → Sidak post-hoc
emmeans(model, ~ Genotype * Diet) %>%
  contrast("pairwise", adjust = "sidak")

# ----------------------------------------------------------------------------------------
# If data doesn't meet normality assumptions and homogeneity of variance: Log-transformation
df_matsuda_index_female_full$matsuda_log <- log(df_matsuda_index_female_full$matsuda_index)

# And then re-check using the same Shapiro-Wilk test and Levene's test
shapiro.test(df_matsuda_index_female_full$matsuda_log)
leveneTest(matsuda_log ~ Genotype * Diet, data = df_matsuda_index_female_full)

model_mat <- aov(matsuda_log ~ Genotype * Diet, data = df_matsuda_index_female_full)
summary(model_mat)

emmeans(model_mat, ~ Genotype | Diet) %>%
  contrast("pairwise", adjust = "sidak")
