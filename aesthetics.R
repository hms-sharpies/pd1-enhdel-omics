# aesthetics

# constants -----

GGPLOT_TEXT_SCALE_FACTOR <- 5/14
TEXT_SIZE <- 6
TITLE_SIZE <- 6
PLOT_TITLE_SIZE <- 7
LINE_WIDTH <- 0.2
LEGEND_SIZE <- unit(3/32, "inch")

# theme -----

figure_layout_theme <- function(text_size = TEXT_SIZE,
                                title_size = TITLE_SIZE,
                                plot_title_size = PLOT_TITLE_SIZE,
                                line_width = LINE_WIDTH,
                                legend_size = LEGEND_SIZE) {
  theme_classic() +
    theme(
      line = element_line(size = line_width),
      text = element_text(size = text_size),
      title = element_text(size = title_size),
      plot.title = element_text(size = plot_title_size, hjust = 0.5),
      axis.ticks = element_line(size = line_width),
      legend.key.size = legend_size,
      legend.key.height = legend_size,
      legend.key.width = legend_size,
      legend.spacing = unit(0, "pt")
    )
}

theme_set(theme_classic() + theme(legend.key.size = LEGEND_SIZE,
                                  legend.key.height = LEGEND_SIZE,
                                  legend.key.width = LEGEND_SIZE))

umap_theme <- function() {
  theme(axis.title.x = element_text(hjust = 0, vjust = 2.5),
        axis.title.y = element_text(hjust = 0, vjust = -1.5),
        axis.ticks = element_blank(),
        axis.text = element_blank())
}

# keys -----

genotype_key <- c(
  "W" = "WT",
  "D" = "EnhDel",
  "K" = "PD-1 KO"
)

cluster_num_key <- c(
  "3" = "dividing",
  "0" = "terminal",
  "2" = "transitory",
  "1" = "effector-like",
  "4" = "progenitor"
)

cluster_ids_key <- c(
  "3" = "dividing",
  "0" = "terminally exhausted",
  "2" = "transitory exhausted",
  "1" = "effector-like",
  "4" = "progenitor exhausted"
)

cluster_abbr_key <- c(
  "progenitor" = "prog.",
  "effector-like" = "eff.like",
  "transitory" = "trans.",
  "terminal" = "term.",
  "dividing" = "div."
)

# palettes -----

palette_genotype <- c(
  W = "#000000",
  D = "#0574D8",
  K = "#CE7474"
)

palette_genotype_named <- c(
  EnhDel = "#0574D8",
  "PD-1 KO" = "#CE7474",
  WT = "#000000"
)

palette_orig.ident <- c(
  D1 = "#0574D8",
  D2 = "#055693",
  K1 = "#CE7474",
  K2 = "#844C4C",
  W1 = "#000000",
  W2 = "#424242"
)

palette_cluster_ids <- c(
  "progenitor exhausted" = "#00B295",
  "effector-like" = "#2B98CA",
  "transitory exhausted" = "#8850B9",
  "terminally exhausted" = "#FF715B",
  "dividing" = "#FFC857"
)

palette_cluster <- c(
  "progenitor" = "#00B295",
  "effector-like" = "#2B98CA",
  "transitory" = "#8850B9",
  "terminal" = "#FF715B",
  "dividing" = "#FFC857"
)

palette_cluster_abbr <- c(
  "prog." = "#00B295",
  "eff.like" = "#2B98CA",
  "trans." = "#8850B9",
  "term." = "#FF715B",
  "div." = "#FFC857"
)

palette_social <- c(
  "#FB3958",
  "#FFC838",
  "#6DC993",
  "#458EFF",
  "#DED1C1",
  "#CCCCFF",
  "#AAD450"
)
