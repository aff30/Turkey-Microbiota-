## Beta diversity - October 26

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
library(microViz)
library(pairwiseAdonis)

# ---- 1) Subset and clean phyloseq ----
ps <- readRDS("ps_object/ps-decontam-filtered-counts.rds")

ps <- ps %>%
  subset_samples(Sample_Type == "Ceca") %>%
  tax_fix()

#principal coordinate analysis

#colors
#colors
trts <- c("Control" = "#593827FF", 
          "1X" = "#A1C2ED",
          "2X" = "#322EA6FF")


## Beta Diversity
ps <- tax_fix(ps)
ps <- phyloseq_validate(ps, remove_undetected = TRUE)


ps %>% 
  tax_transform(trans = "identity", rank = "unique") %>% 
  dist_calc("aitchison") 

sample_data(ps)

# Plot a PCA
# Principal Components Analysis is an unconstrained method that does not use a distance matrix.
# Each point is a sample, and samples that appear closer together are typically 
# more similar to each other than samples which are further apart.

# roughly estimate which samples will contain more of that taxon 
ps %>%
  ps_mutate(Days = factor(as.character(Days),
                          levels = c("38", "49", "97", "126"))) %>%
  tax_fix() %>% 
  tax_transform("clr", rank = "Genus") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>% 
  ord_plot(color = "Days", plot_taxa = 1:7, size = 6, 
           taxa.text.size = 8, parse = TRUE) +
  stat_ellipse(aes(group = Days, color=Days, linetype = Days)) +
  scale_linetype_manual(values = c("dashed", "solid", "dashed", "solid"))+
  scale_color_manual(values = c("#875E3BFF", "#E37F3EFF", "#B3C6ADFF", "#384630FF")) +
  theme_classic() +
  theme(axis.text= element_text(size = 22),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22)
  )


## PCoA
# Aitchison distance
ps %>%
  ps_mutate(Days = factor(as.character(Days),
                          levels = c("38", "49", "97", "126"))) %>%
  tax_fix() %>% 
  tax_transform("clr", rank = "Genus") %>% 
  #tax_transform("identity", rank = "Class") %>% # don't transform!
  #dist_calc("aitchison") %>% 
  ord_calc() %>% 
  ord_plot(color = "Treatment", shape = "Days", size = 7) +
  stat_ellipse(aes(group = Treatment, color=Treatment, linetype = Treatment)) +
  scale_linetype_manual(
    name = "Treatment",
    values = c("Control" = "dashed", "1X" = "solid", "2X" = "dashed"),
    breaks = c("Control", "1X", "2X"),
    labels = c("Control", "CON + 0.0125% DFM", "CON + 0.025% DFM")
  ) +
  scale_color_manual(
    name = "Treatment",
    values = trts,
    breaks = c("Control", "1X", "2X"),
    labels = c("Control", "CON + 0.0125% DFM", "CON + 0.025% DFM")
  ) +
  theme_classic() +
  theme(axis.text= element_text(size = 22),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22)
  )


#PERMANOVA overall trt and age 
# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_fix() %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

#ADONIS test
# age, treatment

vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment, permutations = 10000) 
# p = 0.6223
vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Days, permutations = 10000)
# p = 9.999e-05 ***

# Pairwise comparison using pairwiseAdonis 
data <- phyloseq::sample_data(beta1)
data
pairwiseAdonis::pairwise.adonis(psdist, data$Days)

#pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
#1 126 vs 38  1 1742.0312 14.118131 0.19576397   0.001      0.006   *
#2 126 vs 49  1 1059.6784  8.193496 0.15118965   0.001      0.006   *
#3 126 vs 97  1  414.8707  3.227788 0.06556841   0.001      0.006   *
#4  38 vs 49  1  584.1759  4.628035 0.07389717   0.001      0.006   *
#5  38 vs 97  1 1069.4819  8.515631 0.12802450   0.001      0.006   *
#6  49 vs 97  1  732.2158  5.542614 0.10753459   0.001      0.006   *



# testing if age is masking treatment
# ADONIS with both Age and Treatment
adonis_result <- vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment + phyloseq::sample_data(beta1)$Days, permutations = 10000)
print(adonis_result)

#Df SumOfSqs      R2      F    Pr(>F)    
#Model     5   2085.9 0.19153 4.3117 9.999e-05 ***
# Residual 91   8804.9 0.80847                     
#Total    96  10890.9 1.00000  

vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment * phyloseq::sample_data(beta1)$Days, permutations = 10000)
#Df SumOfSqs      R2      F    Pr(>F)    
#Model    11   2676.0 0.24571 2.5172 9.999e-05 ***
# Residual 85   8214.9 0.75429                     
#Total    96  10890.9 1.00000

#The significant p-value suggests that the interaction between Treatment and Days is meaningful in explaining variation in the mycobiome data.

# Assuming 'psdist' is your distance matrix and the groups are based on Age_Days and Treatment
betadisper_result <- betadisper(psdist, group = interaction(phyloseq::sample_data(beta1)$Days, phyloseq::sample_data(beta1)$Treatment))

# To view the centroids and distances:
betadisper_result$centroids  # this gives the centroid locations for each group
betadisper_result$distances  # this gives the distances to the centroid for each sample

# Assuming 'betadisper_result$centroids' is a matrix or data frame with centroid coordinates
centroids <- betadisper_result$centroids
dist_centroids <- dist(centroids)
dist_centroids

###Beta for the Litter samples###

# ---- 1) Subset and clean phyloseq ----
ps <- readRDS("ps_object/ps-decontam-filtered-counts.rds")

ps <- ps %>%
  subset_samples(Sample_Type == "Litter") %>%
  subset_samples(Days == "34" | Days == "102") %>%
  subset_samples(Treatment == "Control" | Treatment == "1X" | Treatment == "2X") %>%
  tax_fix()

sample_data(ps)


#principal coordinate analysis

#colors
#colors
trts <- c("Control" = "#593827FF", 
          "1X" = "#A1C2ED",
          "2X" = "#322EA6FF")


## Beta Diversity
ps <- tax_fix(ps)
ps <- phyloseq_validate(ps, remove_undetected = TRUE)


ps %>% 
  tax_transform(trans = "identity", rank = "unique") %>% 
  dist_calc("aitchison") 

sample_data(ps)

# Plot a PCA
# Principal Components Analysis is an unconstrained method that does not use a distance matrix.
# Each point is a sample, and samples that appear closer together are typically 
# more similar to each other than samples which are further apart.

# roughly estimate which samples will contain more of that taxon 
ps %>%
  ps_mutate(Days = factor(as.character(Days),
                          levels = c("34", "102"))) %>%
  
  tax_fix() %>% 
  tax_transform("clr", rank = "Genus") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>% 
  ord_plot(color = "Days", plot_taxa = 1:7, size = 8, 
           taxa.text.size = 10, parse = TRUE) +
  stat_ellipse(aes(group = Days, color=Days, linetype = Days)) +
  scale_linetype_manual(values = c("dashed", "solid", "dashed", "solid"))+
  scale_color_manual(values = c("#C33C54", "#254E70")) +
  theme_classic() +
  theme(axis.text= element_text(size = 22),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22)
  )


## PCoA
# Aitchison distance
ps %>%
  ps_mutate(Days = factor(as.character(Days),
                          levels = c("34", "102"))) %>%
  tax_fix() %>% 
  tax_transform("clr", rank = "Genus") %>% 
 # tax_transform("identity", rank = "Class") %>% # don't transform!
  #dist_calc("aitchison") %>% 
  ord_calc() %>% 
  ord_plot(color = "Treatment", shape = "Days", size = 7) +
  stat_ellipse(aes(group = Treatment, color=Treatment, linetype = Treatment)) +
  scale_linetype_manual(
    name = "Treatment",
    values = c("Control" = "dashed", "1X" = "solid", "2X" = "dashed"),
    breaks = c("Control", "1X", "2X"),
    labels = c("Control", "CON + 0.0125% DFM", "CON + 0.025% DFM")
  ) +
  scale_color_manual(
    name = "Treatment",
    values = trts,
    breaks = c("Control", "1X", "2X"),
    labels = c("Control", "CON + 0.0125% DFM", "CON + 0.025% DFM")
  ) +
  theme_classic() +
  theme(axis.text= element_text(size = 22),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22)
  )



#PERMANOVA overall trt and age 
# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_fix() %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

#ADONIS test
# age, treatment

vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment, permutations = 10000) 
# p = 0.99
vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Days, permutations = 10000)
# p = 9.999e-05 ***

# Pairwise comparison using pairwiseAdonis 
data <- phyloseq::sample_data(beta1)
data
pairwiseAdonis::pairwise.adonis(psdist, data$Days)

#pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
#34 vs 102  1  782.5783 6.065707 0.09467947   0.001      0.001  **


# testing if age is masking treatment
# ADONIS with both Age and Treatment
adonis_result <- vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment + phyloseq::sample_data(beta1)$Days, permutations = 10000)
print(adonis_result)


vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment * phyloseq::sample_data(beta1)$Days, permutations = 10000)
#Df SumOfSqs      R2      F Pr(>F)   
#Model     5   1145.3 0.13857 1.7373 0.0024 **
 # Residual 54   7120.2 0.86143                 
#Total    59   8265.6 1.00000                 

#The significant p-value suggests that the interaction between Treatment and Days is meaningful in explaining variation in the mycobiome data.

# Assuming 'psdist' is your distance matrix and the groups are based on Age_Days and Treatment
betadisper_result <- betadisper(psdist, group = interaction(phyloseq::sample_data(beta1)$Days, phyloseq::sample_data(beta1)$Treatment))

# To view the centroids and distances:
betadisper_result$centroids  # this gives the centroid locations for each group
betadisper_result$distances  # this gives the distances to the centroid for each sample

# Assuming 'betadisper_result$centroids' is a matrix or data frame with centroid coordinates
centroids <- betadisper_result$centroids
dist_centroids <- dist(centroids)
dist_centroids
