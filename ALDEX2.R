#Differential relative abundance

set.seed(200789)
require(tidyverse)
require(phyloseq)
library(ALDEx2)
library(microViz)
library(writexl)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# load data
ps <- readRDS("ps_object/ps-decontam-filtered-counts.rds")
ps


ps <- subset_samples(ps, Treatment == "Control" | Treatment == "2X" | Treatment == "1X")
ps <- subset_samples(ps, Sample_Type == "Ceca")

sample_data(ps)


# Inspect tables
head(otu_table(ps), n=3)
head(tax_table(ps), n=3)
head(sample_data(ps), n=3)

#transform the ps 
otu_table <- otu_table(ps)
#transposed_otu_table <- t(otu_table) # check if it need to be transposed 
#otu_table(ps) <- transposed_otu_table

dim(otu_table(ps))  # Check dimensions of the OTU table
head(otu_table(ps))  # View the first few rows of the OTU table

sample_data(ps)

#Two treatments only - it is a t-test  - 
##select two samples
ps <- subset_samples(ps, Days == "38" | Days == "126")

sample_data(ps)

# aggregate at genus level for counts
ps1 <- ps %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

head(sample_data(ps1))

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(ps)$Days, taxa_rank = "all", gamma = 0.5, method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 1000, denom = "none")

sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)
## 

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#Differences between 126 and 38

# make a table of significant corrected p-values
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
sig_aldex2 <- left_join(sig_aldex2, taxa_info)
save(sig_aldex2, file = "tables/table_sig_taxa_age_38_126.RData")
write.table(sig_aldex2, file = "table_sig_taxa_age_38_126.csv", sep = ",", col.names = TRUE, row.names = FALSE)
write_xlsx(sig_aldex2, "tables/table_sig_taxa_age_38_126_TRT_1x.xlsx")


##graph 1
dat1 <- readxl::read_excel("tables/table_sig_taxa_age_38_126_TRT_1x.xlsx")

ggplot(data=dat1, aes(x=Genus, y=(effect), col=Days, label=Genus)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  #scale_color_brewer() +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=0, col="black")


##cleavand dot plot graph
dat1$Days <- factor(dat1$Days)

p1 <- ggplot(dat1, aes(effect, Genus)) +
  geom_line(aes(group = Genus)) +
  geom_point(aes(color = factor(Days)), size = 6) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "gray30") +  # add 0.0 line
  scale_color_manual(values = c("#875E3BFF","#384630FF")) +
  scale_shape_manual(values = c(16, 17)) +
  labs(
    x = "Effect",
    y = "Genus",
    color = "Age (Days)"
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid = element_blank(),          # remove gridlines
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18, face = "italic"),
    legend.title = element_text(size = 20)
  )

p2 <- ggpar(p1, font.ytickslab = "italic")
p2


#Differential relative abundance - Litter

set.seed(200789)
require(tidyverse)
require(phyloseq)
library(ALDEx2)
library(microViz)
library(writexl)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# load data
ps <- readRDS("ps_object/ps-decontam-filtered-counts.rds")
ps


ps <- subset_samples(ps, Treatment == "Control" | Treatment == "2X" | Treatment == "1X")
ps <- subset_samples(ps, Sample_Type == "Litter")

sample_data(ps)


# Inspect tables
head(otu_table(ps), n=3)
head(tax_table(ps), n=3)
head(sample_data(ps), n=3)

#transform the ps 
otu_table <- otu_table(ps)
#transposed_otu_table <- t(otu_table) # check if it need to be transposed 
#otu_table(ps) <- transposed_otu_table

dim(otu_table(ps))  # Check dimensions of the OTU table
head(otu_table(ps))  # View the first few rows of the OTU table

sample_data(ps)

#Two treatments only - it is a t-test  - 
##select two samples
ps <- subset_samples(ps, Days == "34" | Days == "102")

sample_data(ps)

# aggregate at genus level for counts
ps1 <- ps %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

head(sample_data(ps1))

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(ps)$Days, taxa_rank = "all", gamma = 0.5, method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 1000, denom = "none")

sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)
#No differences 


