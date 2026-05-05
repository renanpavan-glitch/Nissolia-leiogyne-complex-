# ==============================================================================
# SCRIPT: NORMALITY TESTS, INFERENTIAL STATISTICS, AND CORRELATIONS
# ==============================================================================

library(Hmisc)
library(FSA)
library(readr)
library(ggplot2)
library(stats)
library(corrplot)
library(tidyverse)
library(ggstatsplot)
library(reshape2)
options(stringsAsFactors = FALSE)

# -------------------- Data Reading and Processing -------------------- #
###---Replace with leaf_fruit.CSV for leaf and fruit analysis---###
famd <- read_csv2("flower.CSV", locale = locale(encoding = "UTF-8")) |> type_convert()

famd[] <- lapply(famd, function(x) {
  if (is.character(x)) {
    as.factor(x)
  } else if (is.numeric(x) && length(unique(x)) <= 3) {
    as.factor(x)
  } else {
    x
  }
})

# ---------------------- Quantitative Data Tests ---------------------- #
famd_numeric <- famd[ , sapply(famd, is.numeric)]

# Normality test (Shapiro-Wilk)
shapiro_results <- apply(famd_numeric, 2, function(x) {
  x <- na.omit(x)
  if (length(x) < 3) return(NA)
  shapiro.test(x)$p.value
})

# Separate variables
normal_vars <- names(shapiro_results[shapiro_results > 0.05])
non_normal_vars <- names(shapiro_results[shapiro_results <= 0.05])

results <- list()

# ANOVA loop for normal variables
for (var in normal_vars) {
  formula <- as.formula(paste0("`", var, "` ~ `Sp.`"))  
  model <- aov(formula, data = famd)
  p <- summary(model)[[1]][["Pr(>F)"]][1]
  posthoc <- TukeyHSD(model)
  
  results[[var]] <- list(
    type = "ANOVA",
    p_value = p,
    posthoc = posthoc
  )
}

# Kruskal-Wallis for non-normal variables
for (var in non_normal_vars) {
  formula <- as.formula(paste0("`", var, "` ~ `Sp.`"))
  kruskal <- kruskal.test(formula, data = famd)
  posthoc <- tryCatch(
    dunnTest(formula, data = famd, method = "bonferroni"),
    error = function(e) NULL
  )
  
  results[[var]] <- list(
    type = "Kruskal-Wallis",
    p_value = kruskal$p.value,
    posthoc = posthoc
  )
}

# SUMMARY
for (var in names(results)) {
  cat("\n=============================\n")
  cat("Variable:", var, "\n")
  cat("Test type:", results[[var]]$type, "\n")
  cat("Overall p-value:", results[[var]]$p_value, "\n")
}

# -------- MIXED CORRELATION MATRIX (Pearson/Spearman) -------- #
n <- ncol(famd_numeric)
cor_matrix <- matrix(NA, n, n)
p_matrix <- matrix(NA, n, n)
type_matrix <- matrix(NA, n, n)
colnames(cor_matrix) <- rownames(cor_matrix) <- colnames(famd_numeric)
colnames(p_matrix) <- rownames(p_matrix) <- colnames(famd_numeric)
colnames(type_matrix) <- rownames(type_matrix) <- colnames(famd_numeric)

for (i in 1:n) {
  for (j in 1:n) {
    x <- famd_numeric[[i]]
    y <- famd_numeric[[j]]
    
    normal_i <- colnames(famd_numeric)[i] %in% normal_vars
    normal_j <- colnames(famd_numeric)[j] %in% normal_vars
    
    if (normal_i && normal_j) {
      test <- cor.test(x, y, method = "pearson")
      type_matrix[i, j] <- "P"
    } else {
      test <- cor.test(x, y, method = "spearman", exact = FALSE)
      type_matrix[i, j] <- "S"
    }
    
    cor_matrix[i, j] <- test$estimate
    p_matrix[i, j] <- test$p.value
  }
}

df_cor <- melt(cor_matrix)
df_type <- melt(type_matrix)
df_p <- melt(p_matrix)

df_all <- cbind(df_cor, type = df_type$value, pval = df_p$value)
colnames(df_all) <- c("Var1", "Var2", "cor", "type", "pval")
df_all$Var1 <- as.character(df_all$Var1)
df_all$Var2 <- as.character(df_all$Var2)
df_all_filtered <- df_all[df_all$Var1 >= df_all$Var2, ]

ggplot(df_all_filtered, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#E46726", mid = "white", high = "#6D9EC1", midpoint = 0,
                       name = "Correlation") +
  geom_text(aes(label = round(cor, 2)), size = 5) +
  geom_text(aes(label = type), nudge_y = -0.25, size = 3, fontface = "italic", color = "black") +
  theme_minimal(base_family = "serif") +  
  labs(title = "Correlation Matrix with Type (P = Pearson, S = Spearman)",
       x = "", y = "") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10))


# ---------------------- Qualitative Data Tests ---------------------- #
names(famd)[names(famd) == "Sp."] <- "Taxa"

morphological_vars <- c("Infl.", "Std._Shape", "Indum.")

for (var in morphological_vars) {
  print(
    ggbarstats(
      data = famd,
      x = Taxa,                 
      y = !!sym(var),            
      title = paste("Association between", var, "and Taxa"),
      results.subtitle = TRUE
    )
  )
}