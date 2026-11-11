##Relative abundance 
#Ana Fonseca 
#September 2026

library(phyloseq)
library(microViz)
library(dplyr)
library(forcats)
library(scales)
library(ggplot2)

# ---- 1) Subset and clean phyloseq ----
ps <- readRDS("ps_object/ps-decontam-filtered-counts.rds")

ps <- ps %>%
  subset_samples(Sample_Type == "Ceca") %>%
  # subset_samples(Days %in% c("38","49")) %>% 
  tax_fix()

# ---- 2) Define hue/shade ranks ----
hueRank       <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank     <- "Genus"  # change to "Genus" if you prefer

# ---- 3) Sort and aggregate at the LOWER rank (Family) ----
pseq2 <- ps %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# ---- 4) Build hierarchical palette info ----
nHues   <- 3  # number of Phyla with distinct hues
nShades <- 4  # number of Families per Phylum with distinct shades

hierarchicalPalInfo <- data.frame(
  hue    = as.vector(tt_get(pseq2)[, hueRank]),
  shade  = as.vector(tt_get(pseq2)[, shadeRank]),
  counts = taxa_sums(otu_get(pseq2))
) %>%
  mutate(
    hue = forcats::fct_other(
      f = hue,
      keep = unique(hue)[seq_len(min(nHues, length(unique(hue))))],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)),
    padHue  = max(nChrHue) - nChrHue
  ) %>%
  group_by(hue) %>%
  mutate(
    shade = forcats::fct_other(
      f = shade,
      keep = unique(shade)[seq_len(min(nShades - 1, length(unique(shade))))],
      other_level = "Other"
    )
  ) %>%
  ungroup() %>%
  mutate(
    nChrShade = nchar(as.character(shade)),
    padShade  = max(nChrShade) - nChrShade,
    Taxa      = paste0(
      hue, ": ", strrep(" ", padHue),
      shade, strrep(" ", padShade)
    )
  )

# unique combined labels for Phylum:Family
taxaNames <- unique(hierarchicalPalInfo$Taxa)

# ---- 5) Make colour matrix and properly NAMED palette ----
hierarchicalPalMatrix <- matrix(
  data = sapply(
    X = seq(from = 30, to = 75, length.out = nShades),
    FUN = function(l) scales::hue_pal(l = l, h.start = 120)(n = nHues)
  ),
  byrow = TRUE,
  ncol  = nHues
)

hierarchicalPalMatrix <- cbind(
  hierarchicalPalMatrix,
  grey.colors(n = nShades)
)

# Flatten, then KEEP only as many colours as taxa
flatCols <- as.vector(hierarchicalPalMatrix)
flatCols <- flatCols[seq_along(taxaNames)]

hierarchicalPal <- setNames(flatCols, taxaNames)

# Now this should work:
# tax_palette_plot(hierarchicalPal) +
#   theme(axis.text.y.left = element_text(family = "mono"))

#Just ceca
# ---- 6) Plot with comp_barplot, turning OFF seriation ----
pseq2 %>%
  ps_get() %>%
  tax_mutate(`Phylum: Genus` = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level    = "Phylum: Genus",
    n_taxa       = length(hierarchicalPal),
    tax_order    = "asis",   # no hclust/bray ordering
    sample_order = "asis",   # <- IMPORTANT: avoid vegdist / seriation
    palette      = hierarchicalPal,
    bar_width    = 0.975
  ) +
  coord_flip() +
  theme(
    legend.text = element_text(family = "mono"),
    axis.text.y  = element_blank(), # remove y labels (sample IDs)
    axis.ticks.y = element_blank()
  )


# ---- 6) Prepare the data once, keeping the same taxa and colours ----
plot_ps <- pseq2 %>%
  ps_get() %>%
  ps_mutate(
    Treatment = dplyr::recode(
      as.character(Treatment),
      "Control" = "Control",
      "1X"      = "CON + 0.0125%",
      "2X"      = "CON + 0.025%"
    ),
    Days = factor(as.character(Days), levels = c("38", "49", "97", "126"))
  ) %>%
  tax_mutate(
    `Phylum: Genus` = hierarchicalPalInfo$Taxa,
    .keep = "none"
  )

# ---- 7) Function: one treatment, days arranged horizontally ----
plot_treatment <- function(treatment_name) {
  plot_ps %>%
    ps_filter(Treatment == treatment_name) %>%
    comp_barplot(
      tax_level    = "Phylum: Genus",
      n_taxa       = length(hierarchicalPal),
      tax_order    = "asis",
      sample_order = "asis",
      palette      = hierarchicalPal,
      bar_width    = 1
    ) +
    coord_flip() +
    facet_grid(
      rows = vars(Days),
      scales = "free_y",
      space = "free_y"
    ) +
    labs(
      title = treatment_name,
      x = NULL,
      y = "Relative abundance",
      fill = "Phylum: Genus"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text = element_text(size = 18, face = "bold"),
      panel.spacing.x = grid::unit(10, "pt"),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.ticks.y = element_blank(),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 15, face = "bold")
    )
}

p_control <- plot_treatment("Control")
p_1x      <- plot_treatment("CON + 0.0125%")
p_2x      <- plot_treatment("CON + 0.025%")

p_control
p_1x
p_2x


library(ggpubr)

combined_plot <- ggarrange(
  p_control, p_1x, p_2x,
  ncol = 3,
  nrow = 1,
  labels = NULL,
  common.legend = TRUE,
  legend = "bottom"
)

combined_plot

p_control_no_days <- p_control +
  theme(
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )

p_1x_no_days <- p_1x +
  theme(
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )

combined_plot <- ggpubr::ggarrange(
  p_control_no_days, p_1x_no_days, p_2x,
  ncol = 3,
  common.legend = TRUE,
  legend = "bottom"
)

combined_plot

legend_settings <- list(
  guides(fill = guide_legend(ncol = 3, byrow = TRUE)),
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.size = grid::unit(0.35, "cm"),
    legend.spacing.x = grid::unit(0.25, "cm")
  )
)

p_control_final <- p_control +
  legend_settings +
  theme(
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )

p_1x_final <- p_1x +
  legend_settings +
  theme(
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )

p_2x_final <- p_2x + legend_settings

combined_plot <- ggpubr::ggarrange(
  p_control_final, p_1x_final, p_2x_final,
  ncol = 3,
  common.legend = TRUE,
  legend = "bottom"
)

combined_plot


###NOW LITTER


library(phyloseq)
library(microViz)
library(dplyr)
library(forcats)
library(scales)
library(ggplot2)

# ---- 1) Subset and clean phyloseq ----
ps <- readRDS("ps_object/ps-decontam-filtered-counts.rds")

ps <- ps %>%
  subset_samples(Sample_Type == "Litter") %>%
  subset_samples(Days %in% c("34","102")) %>% 
  tax_fix()

ps <- prune_samples(sample_names(ps) != "L21", ps)
sample_data(ps)

# ---- 2) Define hue/shade ranks ----
hueRank       <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank     <- "Genus"  # change to "Genus" if you prefer

# ---- 3) Sort and aggregate at the LOWER rank (Family) ----
pseq2 <- ps %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# ---- 4) Build hierarchical palette info ----
nHues   <- 3  # number of Phyla with distinct hues
nShades <- 4  # number of Families per Phylum with distinct shades

hierarchicalPalInfo <- data.frame(
  hue    = as.vector(tt_get(pseq2)[, hueRank]),
  shade  = as.vector(tt_get(pseq2)[, shadeRank]),
  counts = taxa_sums(otu_get(pseq2))
) %>%
  mutate(
    hue = forcats::fct_other(
      f = hue,
      keep = unique(hue)[seq_len(min(nHues, length(unique(hue))))],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)),
    padHue  = max(nChrHue) - nChrHue
  ) %>%
  group_by(hue) %>%
  mutate(
    shade = forcats::fct_other(
      f = shade,
      keep = unique(shade)[seq_len(min(nShades - 1, length(unique(shade))))],
      other_level = "Other"
    )
  ) %>%
  ungroup() %>%
  mutate(
    nChrShade = nchar(as.character(shade)),
    padShade  = max(nChrShade) - nChrShade,
    Taxa      = paste0(
      hue, ": ", strrep(" ", padHue),
      shade, strrep(" ", padShade)
    )
  )

# unique combined labels for Phylum:Family
taxaNames <- unique(hierarchicalPalInfo$Taxa)

# ---- 5) Make colour matrix and properly NAMED palette ----
hierarchicalPalMatrix <- matrix(
  data = sapply(
    X = seq(from = 30, to = 75, length.out = nShades),
    FUN = function(l) scales::hue_pal(l = l, h.start = 30)(n = nHues)
  ),
  byrow = TRUE,
  ncol  = nHues
)

hierarchicalPalMatrix <- cbind(
  hierarchicalPalMatrix,
  grey.colors(n = nShades)
)

# Flatten, then KEEP only as many colours as taxa
flatCols <- as.vector(hierarchicalPalMatrix)
flatCols <- flatCols[seq_along(taxaNames)]

hierarchicalPal <- setNames(flatCols, taxaNames)

# Now this should work:
tax_palette_plot(hierarchicalPal) +
 theme(axis.text.y.left = element_text(family = "mono"))

#Just Litter - whole data
# ---- 6) Plot with comp_barplot, turning OFF seriation ----
pseq2 %>%
  ps_get() %>%
  tax_mutate(`Phylum: Genus` = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level    = "Phylum: Genus",
    n_taxa       = length(hierarchicalPal),
    tax_order    = "asis",   # no hclust/bray ordering
    sample_order = "asis",   # <- IMPORTANT: avoid vegdist / seriation
    palette      = hierarchicalPal,
    bar_width    = 0.975
  ) +
  coord_flip() +
  theme(
    legend.text = element_text(family = "mono"),
    axis.text.y  = element_blank(), # remove y labels (sample IDs)
    axis.ticks.y = element_blank()
  )

sample_data(pseq2)

# ---- 6) Prepare the data once, keeping the same taxa and colours ----
plot_ps <- pseq2 %>%
  ps_get() %>%
  ps_mutate(
    Treatment = dplyr::recode(
      as.character(Treatment),
      "Control" = "Control",
      "1X"      = "CON + 0.0125%",
      "2X"      = "CON + 0.025%"
    ),
    Days = factor(as.character(Days), levels = c("34", "102"))
  ) %>%
  tax_mutate(
    `Phylum: Genus` = hierarchicalPalInfo$Taxa,
    .keep = "none"
  )

# ---- 7) Function: one treatment, days arranged horizontally ----
plot_treatment <- function(treatment_name) {
  plot_ps %>%
    ps_filter(Treatment == treatment_name) %>%
    comp_barplot(
      tax_level    = "Phylum: Genus",
      n_taxa       = length(hierarchicalPal),
      tax_order    = "asis",
      sample_order = "asis",
      palette      = hierarchicalPal,
      bar_width    = 1
    ) +
    coord_flip() +
    facet_grid(
      rows = vars(Days),
      scales = "free_y",
      space = "free_y"
    ) +
    labs(
      title = treatment_name,
      x = NULL,
      y = "Relative abundance",
      fill = "Phylum: Genus"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text = element_text(size = 18, face = "bold"),
      panel.spacing.x = grid::unit(10, "pt"),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.ticks.y = element_blank(),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 15, face = "bold")
    )
}

p_control <- plot_treatment("Control")
p_1x      <- plot_treatment("CON + 0.0125%")
p_2x      <- plot_treatment("CON + 0.025%")

p_control
p_1x
p_2x


library(ggpubr)

combined_plot <- ggarrange(
  p_control, p_1x, p_2x,
  ncol = 3,
  nrow = 1,
  labels = NULL,
  common.legend = TRUE,
  legend = "bottom"
)

combined_plot

p_control_no_days <- p_control +
  theme(
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )

p_1x_no_days <- p_1x +
  theme(
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )

combined_plot <- ggpubr::ggarrange(
  p_control_no_days, p_1x_no_days, p_2x,
  ncol = 3,
  common.legend = TRUE,
  legend = "bottom"
)

combined_plot

legend_settings <- list(
  guides(fill = guide_legend(ncol = 3, byrow = TRUE)),
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.size = grid::unit(0.35, "cm"),
    legend.spacing.x = grid::unit(0.25, "cm")
  )
)

p_control_final <- p_control +
  legend_settings +
  theme(
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )

p_1x_final <- p_1x +
  legend_settings +
  theme(
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )

p_2x_final <- p_2x + legend_settings

combined_plot <- ggpubr::ggarrange(
  p_control_final, p_1x_final, p_2x_final,
  ncol = 3,
  common.legend = TRUE,
  legend = "bottom"
)

combined_plot
