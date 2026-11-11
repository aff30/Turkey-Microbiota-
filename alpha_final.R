#Alpha diversity 
#Ana fonseca 
#September 2026

#load packages 
#generally useful ones 
library(BiocManager)
library(stringr)
library(tidyverse)
library(tidylog)

#ones we specifically need 
library(phyloseq)
library(ggplot2)
library(vegan) # version 2.6-4
library(decontam)
library(microViz)
library(ggpubr)
library(car)

#Load data
ps <- readRDS("ps_object/ps-decontam-filtered-counts.rds")
ps

sample_data(ps)

ps <- subset_samples(ps, Treatment == "Control" | Treatment == "1X" | Treatment == "2X")
ps <- subset_samples(ps, Sample_Type == "Ceca")



# get sample data
sampdf <- sample_data(ps) %>% 
  data.frame() %>% 
  rownames_to_column(var = "Sample_ID") 

# get diversity
alpha <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Observed")) %>% 
  # make ID
  rownames_to_column(var = "Sample_ID") 
#alpha$Sample_ID <- str_remove(alpha$Sample_ID, "X")
#str_remove(alpha$Project_ID,"X")
# merge(sampdf, alpha, by = "Project_ID") 
alpha <- merge(sampdf, alpha, by = "Sample_ID") 

## ---- Shannon Kruskal Wallis  ----
#test for normality 
ggdensity(alpha$Shannon, 
          main = "Density plot of diversity",
          xlab = "diversity")

ggqqplot(alpha$Shannon)
qqPlot(alpha$Shannon)

#Not normal but just to confirm 
shapiro.test(alpha$Shannon)


kruskalShannonTRT <- kruskal.test(Shannon ~ Treatment, data = alpha) 
kruskalShannonTRT
#  p-value = 0.207

kruskalObservedTRT <- kruskal.test(Observed ~ Treatment, data = alpha)
kruskalObservedTRT
# p-value = 0.2477


#Age 
kruskalShannonAGE <- kruskal.test(Shannon ~ Days, data = alpha) 
kruskalShannonAGE
# p-value = 2.9e-05

kruskalObservedAGE <- kruskal.test(Observed ~ Days, data = alpha)
kruskalObservedAGE
# p-value = 4.506e-07


#colors
trts <- c("Control" = "#593827FF", 
          "1X" = "#A1C2ED",
          "2X" = "#322EA6FF")



#boxplots

shannonbox <- 
  ggplot(alpha,
         aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = trts) +
  xlab("") + geom_point(position = "jitter", color = "black", aes(shape = Treatment), size = 3) +
  theme_classic() + 
  labs(title = "",
       x = "", 
       y = "Shannon's Index",
       shape = "Treatment",
       fill = "Treatment") + 
  #scale_x_discrete(name="", labels = c("cloaca" = "Cloaca", "sock"="Sock")) +
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none")

shannonbox

observedbox <- 
  ggplot(alpha, aes(x = Treatment, y = Observed, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = trts) +
  xlab("") +
  geom_point(position = "jitter", color = "black", aes(shape = Treatment), size = 3) +
  theme_classic() + 
  labs(title = "",
       x = "", 
       y = "Observed ASVs",
       shape = "Treatment",
       fill = "Treatment") + 
  #scale_x_discrete(name="", labels = c("cloaca" = "Cloaca", "sock"="Sock")) +
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none")

observedbox


##line plot 
alpha2 <- alpha %>%
  dplyr::mutate(
    Days = factor(as.character(Days), levels = c("38", "49", "97", "126")),
    Treatment = dplyr::recode(
      Treatment,
      "1X" = "CON + 0.0125% DFM",
      "2X" = "CON + 0.025% DFM"
    )
  )


alpha2 %>%
  ggline(
    x = "Days",
    y = "Shannon",
    group = "Treatment",
    color = "Treatment",
    palette = c("#593827FF", "#A1C2ED", "#322EA6FF")
    ,
    add = "mean_sd",
    error.plot = "errorbar",
    linewidth = 1.5,
    add.params = list(size = 3)
  ) +
  theme_classic() +
  font("xy.text",   size = 18) +
  font("xy.title",  size = 18) +
  font("legend.title", size = 16) +
  font("legend.text",  size = 16) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.line    = element_line(linewidth = 1),
    axis.ticks   = element_line(linewidth = 1),
    legend.position = "top"
  )



##Litter alpha
#Load data
ps <- readRDS("ps_object/ps-decontam-filtered-counts.rds")
ps

sample_data(ps)

ps <- subset_samples(ps, Treatment == "Control" | Treatment == "1X" | Treatment == "2X")
ps <- subset_samples(ps, Sample_Type == "Litter")
ps <- subset_samples(ps, Days == "34" | Days == "102")


# get sample data
sampdf <- sample_data(ps) %>% 
  data.frame() %>% 
  rownames_to_column(var = "Sample_ID") 

# get diversity
alpha <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Observed")) %>% 
  # make ID
  rownames_to_column(var = "Sample_ID") 
#alpha$Sample_ID <- str_remove(alpha$Sample_ID, "X")
#str_remove(alpha$Project_ID,"X")
# merge(sampdf, alpha, by = "Project_ID") 
alpha <- merge(sampdf, alpha, by = "Sample_ID") 

## ---- Shannon Kruskal Wallis  ----
#test for normality 
ggdensity(alpha$Shannon, 
          main = "Density plot of diversity",
          xlab = "diversity")

ggqqplot(alpha$Shannon)
qqPlot(alpha$Shannon)

#Not normal but just to confirm 
shapiro.test(alpha$Shannon)


kruskalShannonTRT <- kruskal.test(Shannon ~ Treatment, data = alpha) 
kruskalShannonTRT
#  p-value = 0.6347

kruskalObservedTRT <- kruskal.test(Observed ~ Treatment, data = alpha)
kruskalObservedTRT
# p-value =  0.661

#Age 
kruskalShannonAGE <- kruskal.test(Shannon ~ Days, data = alpha) 
kruskalShannonAGE
# p-value = 0.3815


kruskalObservedAGE <- kruskal.test(Observed ~ Days, data = alpha)
kruskalObservedAGE
# p-value =  0.3119


#colors
trts <- c("Control" = "#593827FF", 
          "1X" = "#A1C2ED",
          "2X" = "#322EA6FF")



#boxplots

shannonbox <- 
  ggplot(alpha,
         aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = trts) +
  xlab("") + geom_point(position = "jitter", color = "black", aes(shape = Treatment), size = 3) +
  theme_classic() + 
  labs(title = "",
       x = "", 
       y = "Shannon's Index",
       shape = "Treatment",
       fill = "Treatment") + 
  #scale_x_discrete(name="", labels = c("cloaca" = "Cloaca", "sock"="Sock")) +
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none")

shannonbox

observedbox <- 
  ggplot(alpha, aes(x = Treatment, y = Observed, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = trts) +
  xlab("") +
  geom_point(position = "jitter", color = "black", aes(shape = Treatment), size = 3) +
  theme_classic() + 
  labs(title = "",
       x = "", 
       y = "Observed ASVs",
       shape = "Treatment",
       fill = "Treatment") + 
  #scale_x_discrete(name="", labels = c("cloaca" = "Cloaca", "sock"="Sock")) +
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none")

observedbox


##line plot 
alpha2 <- alpha %>%
  dplyr::mutate(
    Days = factor(as.character(Days), levels = c("34", "102")),
    Treatment = dplyr::recode(
      Treatment,
      "1X" = "CON + 0.0125% DFM",
      "2X" = "CON + 0.025% DFM"
    )
  )


alpha2 %>%
  ggline(
    x = "Days",
    y = "Shannon",
    group = "Treatment",
    color = "Treatment",
    palette = c("#A1C2ED", "#322EA6FF", "#593827FF")
    ,
    add = "mean_sd",
    error.plot = "errorbar",
    linewidth = 1.5,
    add.params = list(size = 3)
  ) +
  theme_classic() +
  font("xy.text",   size = 18) +
  font("xy.title",  size = 18) +
  font("legend.title", size = 16) +
  font("legend.text",  size = 16) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.line    = element_line(linewidth = 1),
    axis.ticks   = element_line(linewidth = 1),
    legend.position = "top"
  )

