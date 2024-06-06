rm(list = ls())     # clear objects  
graphics.off() 
#######################################
###### Contaminated Burgers ###########
#######################################


# Packages ----------------------------------------------------------------
inst <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("tidyverse","cluster", "factoextra","NbClust","tidyr", 
              "ggplot2", "ggpubr", "broom", "AICcmodavg", "ggcorrplot", 
              "fpc","cluster", "readxl", "magrittr","hrbrthemes",
              "multipanelfigure","klaR","psych","MASS","ggord","devtools",
              "reshape2","RColorBrewer","SensoMineR","FactoMineR","stats",
              "dplyr","writexl","gtools","ggbiplot","ggrepel",
              "ggstatsplot", "plotly", "car", "ez", "openxlsx","reticulate",
              "rstatix")
inst(packages)
theme_set(theme_minimal())


# -------------------------- #
#    Importing Data set      #
# -------------------------- #
(df <- read_excel("Burger.xlsx"))


# -------------------------- #
#    Data exploring          #
# -------------------------- #
dim(df)
str(df)
df <- df %>%
  mutate(
    mo = as.factor(mo),
    time = as.factor(time),
    treatment = as.factor(treatment),
    rep = as.factor(rep),
    Y = as.numeric(Y)
  )
summary(df)


# -------------------------- #
#  Subsets transformations   #
# -------------------------- #
# M.O. subsets
df_A <- df %>% filter(mo == "APC")
df_T <- df %>% filter(mo == "Total_coliforms")
df_E <- df %>% filter(mo == "E_coli")
df_S <- df %>% filter(mo == "S_aureus")


# Negative control mean
CN_A_mean <- mean(df_A$Y[df_A$treatment == "CN"], na.rm = TRUE)
CN_T_mean <- mean(df_T$Y[df_T$treatment == "CN"], na.rm = TRUE)
CN_E_mean <- mean(df_E$Y[df_E$treatment == "CN"], na.rm = TRUE)
CN_S_mean <- mean(df_S$Y[df_S$treatment == "CN"], na.rm = TRUE)


# M.O. subsets with negative control subtracted
df_A <- df_A %>%
  slice(4:n()) %>% 
  mutate(Y =  Y - CN_A_mean)

df_T <- df_T %>%
  slice(4:n()) %>% 
  mutate(Y =  Y - CN_T_mean)

df_E <- df_E %>%
  slice(4:n()) %>% 
  mutate(Y =  Y - CN_E_mean)

df_S <- df_S %>%
  slice(4:n()) %>% 
  mutate(Y =  Y - CN_S_mean)



# -------------------------- #
#    Stat analysis            #
# -------------------------- #
###################
####### APC #######
###################
# dispersion measures
df_A_summary <- df_A %>%
  group_by(time, treatment) %>%
  summarise(
    mean_Y = mean(Y),
    sd_Y = sd(Y),
    n = n(),
    se_Y = sd_Y / sqrt(n),
    .groups = 'drop'
  )


#  ANOVA repeated meassures  #
# -------------------------- #
# First Plot
ggplot(df_A, aes(x = time, y = Y, color = treatment)) +
  geom_boxplot() +
  labs(title = "Distribución de Y por tiempo y tratamiento", x = "Tiempo", y = "UFC/g")


# Plotly Line graph time series
(r_a <- plot_ly(df_A_summary, x = ~time, y = ~mean_Y, color = ~treatment, 
        type = 'scatter', mode = 'lines+markers', 
        error_y = list(type = 'data', array = ~se_Y),
        marker = list(size = 5)) %>%
  layout(title = "Time series APC",
         xaxis = list(title = "Time (days)"),
         yaxis = list(title = "UFC/g"),
         showlegend = TRUE))

htmlwidgets::saveWidget(r_a, "r_a.html")


# Q-Q plot by grups
ggplot(df_A, aes(sample = Y)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(time ~ treatment)

# Shapiro-Wilk test by groups
df_A %>%
  group_by(time, treatment) %>%
  shapiro_test(Y)

# ANOVA and Sphericity check
(modelo_ez_A <- ezANOVA(data = df_A, dv = Y, wid = rep, within = time, between = treatment))


# AOV model Repeated meassures
modelo_aov_A <- aov(Y ~ time * treatment + Error(rep/time), data = df_A)
summary(modelo_aov_A)

# Export data



#   ANOVA time 14 days       #
# -------------------------- #
# Filter
df_A_14 <- df_A %>%
  filter(time == "14") %>% 
  mutate(treatment = factor(treatment, levels = c("C", "PBF", "BU1", "BU2")))

# Perform Shapiro-Wilk test for each treatment group
(shapiro_results_A <- df_A_14 %>%
    group_by(treatment) %>%
    summarise(
      W = shapiro.test(Y)$statistic,
      p_value = shapiro.test(Y)$p.value
    ))

# Perform Levene test for homocedasticity
leveneTest(Y ~ treatment, data = df_A_14)

# Plotting with stats
ggbetweenstats(
  data = df_A_14,
  x = "treatment",
  y = "Y",
  grouping.var = treatment,
  plot.type = "box",
  xlab = "Treatments",
  ylab = "UFC/g",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  title = "APC treatment differences at the end of test",
  type = "parametric", # Normal distribution proved
  var.equal = TRUE # Homocedasticity proved
)

# Extracting Stats data
extract_subtitle(p)
extract_caption(p)
extract_stats(p)

# Join plots
fig <- subplot(
  r_a,
  ggplotly(a_a),  # Transform ggplot2 to plotly
  nrows = 1,      # all plots in 1 row
  shareX = FALSE, # Independent X axis for each plot
  shareY = FALSE  # Independent Y axis for each plot
)
fig



###################
### Coliforms #####
###################
# dispersion measures
df_T_summary <- df_T %>%
  group_by(time, treatment) %>%
  summarise(
    mean_Y = mean(Y),
    sd_Y = sd(Y),
    n = n(),
    se_Y = sd_Y / sqrt(n),
    .groups = 'drop'
  )


#  ANOVA repeated meassures  #
# -------------------------- #
# Plotly Line graph time series
plot_ly(df_T_summary, x = ~time, y = ~mean_Y, color = ~treatment, 
        type = 'scatter', mode = 'lines+markers', 
        error_y = list(type = 'data', array = ~se_Y),
        marker = list(size = 5)) %>%
  layout(title = "Time series Total Coliforms",
         xaxis = list(title = "Time (days)"),
         yaxis = list(title = "UFC/g"),
         showlegend = TRUE)

# Q-Q plot by grups
ggplot(df_T, aes(sample = Y)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(time ~ treatment)

# Shapiro-Wilk test by groups
df_T %>%
  group_by(time, treatment) %>%
  shapiro_test(Y)

# Sphericity check
(modelo_ez_T <- ezANOVA(data = df_T, dv = Y, wid = rep, within = time, between = treatment))

"""# AOV model Repeated meassures alternative
modelo_aov_T <- aov(Y ~ time * treatment + Error(rep/time), data = df_T)
summary(modelo_aov_T)"""


#   ANOVA time 14 days       #
# -------------------------- #
# Filter
df_T_14 <- df_T %>%
  filter(time == "14") %>% 
  mutate(treatment = factor(treatment, levels = c("C", "PBF", "BU1", "BU2")))

# Perform Shapiro-Wilk test for each treatment group
(shapiro_results_T <- df_T_14 %>%
    group_by(treatment) %>%
    summarise(
      W = shapiro.test(Y)$statistic,
      p_value = shapiro.test(Y)$p.value
    ))

# Perform Levene test for homocedasticity
leveneTest(Y ~ treatment, data = df_T_14)

# Plotting with stats
ggbetweenstats(
  data = df_T_14,
  x = "treatment",
  y = "Y",
  grouping.var = treatment,
  plot.type = "box",
  xlab = "Treatments",
  ylab = "UFC/g",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  title = "Total coliforms treatment differences at the end of test",
  type = "parametric", # Normal distribution proved
  var.equal = TRUE # Homocedasticity proved
)

# Extracting Stats data
extract_subtitle(p)
extract_caption(p)
extract_stats(p)



###################
#### E. coli ######
###################
# dispersion measures
df_E_summary <- df_E %>%
  group_by(time, treatment) %>%
  summarise(
    mean_Y = mean(Y),
    sd_Y = sd(Y),
    n = n(),
    se_Y = sd_Y / sqrt(n),
    .groups = 'drop'
  )


#  ANOVA repeated meassures  #
# -------------------------- #
# Plotly Line graph time series
plot_ly(df_E_summary, x = ~time, y = ~mean_Y, color = ~treatment, 
        type = 'scatter', mode = 'lines+markers', 
        error_y = list(type = 'data', array = ~se_Y),
        marker = list(size = 5)) %>%
  layout(title = "Time series <i>E. coli</i>",
         xaxis = list(title = "Time (days)"),
         yaxis = list(title = "UFC/g"),
         showlegend = TRUE)

# Q-Q plot by grups
ggplot(df_E, aes(sample = Y)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(time ~ treatment)

# Shapiro-Wilk test by groups
df_E %>%
  group_by(time, treatment) %>%
  shapiro_test(Y)

# Sphericity check
(modelo_ez_E <- ezANOVA(data = df_E, dv = Y, wid = rep, within = time, between = treatment))


#   ANOVA time 14 days       #
# -------------------------- #
# Filter
df_E_14 <- df_E %>%
  filter(time == "14") %>% 
  mutate(treatment = factor(treatment, levels = c("C", "PBF", "BU1", "BU2")))

# Perform Shapiro-Wilk test for each treatment group
(shapiro_results_E <- df_E_14 %>%
    group_by(treatment) %>%
    summarise(
      W = shapiro.test(Y)$statistic,
      p_value = shapiro.test(Y)$p.value
    ))

# Perform Levene test for homocedasticity
leveneTest(Y ~ treatment, data = df_E_14)

# Plotting with stats
ggbetweenstats(
  data = df_E_14,
  x = "treatment",
  y = "Y",
  grouping.var = treatment,
  plot.type = "box",
  xlab = "Treatments",
  ylab = "UFC/g",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  title = expression(italic("E. coli")~"treatment differences at the end of test"),
  type = "parametric", # Normal distribution proved
  var.equal = TRUE # Homocedasticity proved
)

# Extracting Stats data
extract_subtitle(p)
extract_caption(p)
extract_stats(p)



###################
#### S. aureus ####
###################
# dispersion measures
df_S_summary <- df_S %>%
  group_by(time, treatment) %>%
  summarise(
    mean_Y = mean(Y),
    sd_Y = sd(Y),
    n = n(),
    se_Y = sd_Y / sqrt(n),
    .groups = 'drop'
  )


#  ANOVA repeated meassures  #
# -------------------------- #
# Plotly Line graph time series
plot_ly(df_S_summary, x = ~time, y = ~mean_Y, color = ~treatment, 
        type = 'scatter', mode = 'lines+markers', 
        error_y = list(type = 'data', array = ~se_Y),
        marker = list(size = 5)) %>%
  layout(title = "Time series <i>S. aureus</i>",
         xaxis = list(title = "Time (days)"),
         yaxis = list(title = "UFC/g"),
         showlegend = TRUE)

# Q-Q plot by grups
ggplot(df_S, aes(sample = Y)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(time ~ treatment)

# Shapiro-Wilk test by groups
df_S %>%
  group_by(time, treatment) %>%
  shapiro_test(Y)

# Sphericity check
(modelo_ez_S <- ezANOVA(data = df_S, dv = Y, wid = rep, within = time, between = treatment))


#   ANOVA time 14 days       #
# -------------------------- #
# Filter
df_S_14 <- df_S %>%
  filter(time == "14") %>% 
  mutate(treatment = factor(treatment, levels = c("C", "PBF", "BU1", "BU2")))

# Perform Shapiro-Wilk test for each treatment group
(shapiro_results_S <- df_S_14 %>%
  group_by(treatment) %>%
  summarise(
    W = shapiro.test(Y)$statistic,
    p_value = shapiro.test(Y)$p.value
  ))

# Perform Levene test for homocedasticity
leveneTest(Y ~ treatment, data = df_S_14)

# Plotting with stats
(p <- ggbetweenstats(
  data  = df_S_14,
  x = "treatment",
  y = "Y",
  grouping.var = treatment,
  plot.type = "box",
  xlab = "Treatments",
  ylab = "UFC/g",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  title = expression(italic("S. aureus")~"treatment differences at the end of test"),
  type = "parametric", # Normal distribution proved
  var.equal = TRUE # Homocedasticity proved
))

extract_subtitle(p)
extract_caption(p)
extract_stats(p)