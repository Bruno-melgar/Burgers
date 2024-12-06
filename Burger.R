rm(list = ls())     # clear objects  
graphics.off() 
#ctrl+L   #to clean console
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
              "rstatix", "patchwork")
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
    y = as.numeric(Y)
  ) %>%
  select(-Y)
summary(df)


# -------------------------- #
#  Subsets transformations   #
# -------------------------- #
# M.O. subsets
df_A <- df %>% filter(mo == "APC")
df_T <- df %>% filter(mo == "Total coliforms")
df_E <- df %>% filter(mo == "E. coli")
df_S <- df %>% filter(mo == "S. aureus")
df_M <- df %>% filter(mo == "Salmonella spp")

""" #not used for the second iteration
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
"""


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
    mean_y = mean(y),
    sd_y = sd(y),
    n = n(),
    se_y = sd_y / sqrt(n),
    .groups = 'drop'
  )


#  ANOVA repeated meassures  #
# -------------------------- #
# First Plot
ggplot(df_A, aes(x = time, y = y, color = treatment)) +
  geom_boxplot() +
  labs(title = "Distribución de Y por tiempo y tratamiento", x = "Tiempo", y = "UFC/g")


# Plotly Line graph time series
(fig1.1 <- plot_ly(df_A_summary, x = ~time, y = ~mean_y, color = ~treatment, 
        type = 'scatter', mode = 'lines+markers', 
        error_y = list(type = 'data', array = ~se_y),
        marker = list(size = 5)) %>%
  layout(title = "Time series APC",
         xaxis = list(title = "Time (days)"),
         yaxis = list(title = ~paste("Log<sub>10</sub> CFU/g"), range = c(5, 9.2)),
         showlegend = TRUE))



htmlwidgets::saveWidget(r_a, "r_a.html")


# Q-Q plot by grups
ggplot(df_A, aes(sample = y)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(time ~ treatment)

# Shapiro-Wilk test by groups
df_A %>%
  group_by(time, treatment) %>%
  shapiro_test(y)

# ANOVA and Sphericity check
(modelo_ez_A <- ezANOVA(data = df_A, dv = y, wid = rep, within = time, between = treatment))


# AOV model Repeated meassures
modelo_aov_A <- aov(y ~ time * treatment + Error(rep/time), data = df_A)
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
      W = shapiro.test(y)$statistic,
      p_value = shapiro.test(y)$p.value
    ))

# Perform Levene test for homocedasticity
leveneTest(y ~ treatment, data = df_A_14)

# Plotting with stats
(fig1.2 <- ggbetweenstats(
  data = df_A_14,
  x = "treatment",
  y = "y",
  grouping.var = treatment,
  plot.type = "box",
  xlab = "Treatments",
  ylab = expression(Log[10] ~ "CFU/g"),
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  title = "APC treatment differences at the end of test",
  type = "parametric", # Normal distribution proved
  var.equal = TRUE # Homocedasticity proved
)) 

# Extracting Stats data
extract_subtitle(p)
extract_caption(p)
extract_stats(p)



###################
### Coliforms #####
###################
# dispersion measures
df_T_summary <- df_T %>%
  group_by(time, treatment) %>%
  summarise(
    mean_y = mean(y),
    sd_y = sd(y),
    n = n(),
    se_y = sd_y / sqrt(n),
    .groups = 'drop'
  )


#  ANOVA repeated meassures  #
# -------------------------- #
# Plotly Line graph time series
(fig2.1 <- plot_ly(df_T_summary, x = ~time, y = ~mean_y, color = ~treatment, 
        type = 'scatter', mode = 'lines+markers', 
        error_y = list(type = 'data', array = ~se_y),
        marker = list(size = 5)) %>%
  layout(title = "Time series Total Coliforms",
         xaxis = list(title = "Time (days)"),
         yaxis = list(title = ~paste("Log<sub>10</sub> CFU/g"), range = c(1, 4.2)),
         showlegend = TRUE))

# Q-Q plot by grups
ggplot(df_T, aes(sample = y)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(time ~ treatment)

# Shapiro-Wilk test by groups
df_T %>%
  group_by(time, treatment) %>%
  shapiro_test(y)

df_T_filtered <- df_T %>%
  group_by(time, treatment) %>%
  filter(n_distinct(y) > 1)  # Keep only groups with more than one unique y value

df_T_filtered %>%
  shapiro_test(y) 

# Sphericity check
(modelo_ez_T <- ezANOVA(data = df_T, dv = y, wid = rep, within = time, between = treatment))

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
      W = shapiro.test(y)$statistic,
      p_value = shapiro.test(y)$p.value
    ))

# Perform Levene test for homocedasticity
leveneTest(y ~ treatment, data = df_T_14)

# Plotting with stats
(fig2.2 <- ggbetweenstats(
  data = df_T_14,
  x = "treatment",
  y = "y",
  grouping.var = treatment,
  plot.type = "box",
  xlab = "Treatments",
  ylab = expression(Log[10] ~ "CFU/g"),
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  title = "Total coliforms treatment differences at the end of test",
  type = "parametric", # Normal distribution proved
  var.equal = TRUE # Homocedasticity proved
)) 

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
    mean_y = mean(y),
    sd_y = sd(y),
    n = n(),
    se_y = sd_y / sqrt(n),
    .groups = 'drop'
  )


#  ANOVA repeated meassures  #
# -------------------------- #
# Plotly Line graph time series
(fig3.1 <- plot_ly(df_E_summary, x = ~time, y = ~mean_y, color = ~treatment, 
        type = 'scatter', mode = 'lines+markers', 
        error_y = list(type = 'data', array = ~se_y),
        marker = list(size = 5)) %>%
  layout(title = "Time series <i>E. coli</i>",
         xaxis = list(title = "Time (days)"),
         yaxis = list(title = ~paste("Log<sub>10</sub> CFU/g"), range = c(4, 5.2)),
         showlegend = TRUE))

# Q-Q plot by grups
ggplot(df_E, aes(sample = y)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(time ~ treatment)

# Shapiro-Wilk test by groups
df_E %>%
  group_by(time, treatment) %>%
  shapiro_test(y)

# Sphericity check
(modelo_ez_E <- ezANOVA(data = df_E, dv = y, wid = rep, within = time, between = treatment))


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
      W = shapiro.test(y)$statistic,
      p_value = shapiro.test(y)$p.value
    ))

# Perform Levene test for homocedasticity
leveneTest(y ~ treatment, data = df_E_14)

# Plotting with stats
(fig3.2 <- ggbetweenstats(
  data = df_E_14,
  x = "treatment",
  y = "y",
  grouping.var = treatment,
  plot.type = "box",
  xlab = "Treatments",
  ylab = expression(Log[10] ~ "CFU/g"),
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  title = expression(italic("E. coli")~"treatment differences at the end of test"),
  type = "parametric", # Normal distribution proved
  var.equal = TRUE # Homocedasticity proved
))

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
    mean_y = mean(y),
    sd_y = sd(y),
    n = n(),
    se_y = sd_y / sqrt(n),
    .groups = 'drop'
  )


#  ANOVA repeated meassures  #
# -------------------------- #
# Plotly Line graph time series
(fig4.1 <- plot_ly(df_S_summary, x = ~time, y = ~mean_y, color = ~treatment, 
        type = 'scatter', mode = 'lines+markers', 
        error_y = list(type = 'data', array = ~se_y),
        marker = list(size = 5)) %>%
  layout(title = "Time series <i>S. aureus</i>",
         xaxis = list(title = "Time (days)"),
         yaxis = list(title = ~paste("Log<sub>10</sub> CFU/g"), range = c(4, 5.2)),
         showlegend = TRUE))

# Q-Q plot by grups
ggplot(df_S, aes(sample = y)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(time ~ treatment)

# Shapiro-Wilk test by groups
df_S %>%
  group_by(time, treatment) %>%
  shapiro_test(y)

# Sphericity check
(modelo_ez_S <- ezANOVA(data = df_S, dv = y, wid = rep, within = time, between = treatment))


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
    W = shapiro.test(y)$statistic,
    p_value = shapiro.test(y)$p.value
  ))

# Perform Levene test for homocedasticity
leveneTest(y ~ treatment, data = df_S_14)

# Plotting with stats
(fig4.2 <- ggbetweenstats(
  data  = df_S_14,
  x = "treatment",
  y = "y",
  grouping.var = treatment,
  plot.type = "box",
  xlab = "Treatments",
  ylab = expression(Log[10] ~ "CFU/g"),
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  title = expression(italic("S. aureus")~"treatment differences at the end of test"),
  type = "parametric", # Normal distribution proved
  var.equal = TRUE # Homocedasticity proved
))


plotly::ggplotly(fig4.2, width = 751, height = 508)

extract_subtitle(p)
extract_caption(p)
extract_stats(p)


###################
#### Salmonella ###
###################
# dispersion measures
df_M_summary <- df_M %>%
  group_by(time, treatment) %>%
  summarise(
    mean_y = mean(y),
    sd_y = sd(y),
    n = n(),
    se_y = sd_y / sqrt(n),
    .groups = 'drop'
  )



#  ANOVA repeated meassures  #
# -------------------------- #
# First Plot
ggplot(df_M, aes(x = time, y = y, color = treatment)) +
  geom_boxplot() +
  labs(title = "Distribución de Y por tiempo y tratamiento", x = "Tiempo", y = "UFC/g")


# Plotly Line graph time series
(fig1.1 <- plot_ly(df_M_summary, x = ~time, y = ~mean_y, color = ~treatment, 
                   type = 'scatter', mode = 'lines+markers', 
                   error_y = list(type = 'data', array = ~se_y),
                   marker = list(size = 5)) %>%
    layout(title = "Time series <i>Salmonella</i> spp",
           xaxis = list(title = "Time (days)"),
           yaxis = list(title = ~paste("Log<sub>10</sub> CFU/g"), range = c(1, 6.5)),
           showlegend = TRUE))



htmlwidgets::saveWidget(r_a, "r_a.html")


# Forced interpolation between 7 and 14
df_M_summary <- df_M_summary %>%
  group_by(treatment) %>% # Asegura que las funciones trabajen por grupo
  arrange(time, .by_group = TRUE) %>% # Ordena por tiempo dentro de cada tratamiento
  mutate(mean_y = if_else(time == 10,
                          (lag(mean_y) + lead(mean_y)) / 2, 
                          mean_y),
         sd_y = if_else(time == 10,
                        (lag(sd_y) + lead(sd_y)) / 2, 
                        sd_y)) %>%
  ungroup()


# Plotly Line graph time series corrected by interpolation
(fig1.1 <- plot_ly(df_M_summary, x = ~time, y = ~mean_y, color = ~treatment, 
                   type = 'scatter', mode = 'lines+markers', 
                   error_y = list(type = 'data', array = ~se_y),
                   marker = list(size = 5)) %>%
    layout(title = "Time series <i>Salmonella</i> spp",
           xaxis = list(title = "Time (days)"),
           yaxis = list(title = ~paste("Log<sub>10</sub> CFU/g"), range = c(1, 6.5)),
           showlegend = TRUE))


# Q-Q plot by grups
ggplot(df_M, aes(sample = y)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(time ~ treatment)

# Shapiro-Wilk test by groups
df_M %>%
  group_by(time, treatment) %>%
  shapiro_test(y)

# ANOVA and Sphericity check
(modelo_ez_M <- ezANOVA(data = df_M, dv = y, wid = rep, within = time, between = treatment))


# AOV model Repeated meassures
modelo_aov_M <- aov(y ~ time * treatment + Error(rep/time), data = df_M)
summary(modelo_aov_M)

# Export data



#   ANOVA time 14 days       #
# -------------------------- #
# Filter
df_M_14 <- df_M %>%
  filter(time == "14") %>% 
  mutate(treatment = factor(treatment, levels = c("C", "PBF", "BU1", "BU2")))

# Perform Shapiro-Wilk test for each treatment group
(shapiro_results_M <- df_M_14 %>%
    group_by(treatment) %>%
    summarise(
      W = shapiro.test(y)$statistic,
      p_value = shapiro.test(y)$p.value
    ))

# Perform Levene test for homocedasticity
leveneTest(y ~ treatment, data = df_M_14)

# Plotting with stats
(fig1.2 <- ggbetweenstats(
  data = df_M_14,
  x = "treatment",
  y = "y",
  grouping.var = treatment,
  plot.type = "box",
  xlab = "Treatments",
  ylab = expression(Log[10] ~ "CFU/g"),
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  title = "APC treatment differences at the end of test",
  type = "parametric", # Normal distribution proved
  var.equal = TRUE # Homocedasticity proved
)) 

# Extracting Stats data
extract_subtitle(p)
extract_caption(p)
extract_stats(p)







# -------------------------- #
#    Figures  #
# -------------------------- #
# Mixed panel plot 
# 2 on top, 1 down
# Patchwork package
fig1.2 <- fig1.2 + scale_y_continuous(limits = c(2, 10))
fig2.2 <- fig2.2 + scale_y_continuous(limits = c(2, 10))
fig3.2 <- fig3.2 + scale_y_continuous(limits = c(2, 10))
fig4.2 <- fig4.2 + scale_y_continuous(limits = c(2, 10)) 

combined_plot <- fig1.2 + fig2.2 + fig3.2 + fig4.2 + 
  plot_layout(guides = 'collect')
combined_plot

# Plotly package subplot function
fig <- subplot(fig1.1, fig2.1, fig3.1, fig4.1, nrows = 2, shareX = TRUE, shareY = TRUE) # 2x2 grid
fig <- fig %>% layout(
  title = "Combined Plots", 
  annotations = list(text = "I", x = 0, xref = 'paper', y = 1, yref = 'paper', showarrow = FALSE)) # Annotations
fig



fig1.1 + fig2.1 + fig3.1 + fig4.1
