## ---------------------------------------------------------------------
## Analyzing the source-neuron x downstream-target proportion matrix
## Input: 6-SRC_Proportion_downstream_summary.xlsx  (Proportion_downstream sheet)
##

## ---------------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(pheatmap)
library(grid)
library(circlize)
library(igraph)
library(ggraph)
library(tidygraph)

# ---- 0. Neuron color scheme -------------------------------------------
# One color per neuron/group name, as specified. Both the singular
# ("chordotonal", used as a Source_neuron_name) and plural ("chordotonals",
# used as a downstream-target column name) keys are included so the same
# color applies whichever side of the matrix the name appears on.
neuron_color_map <- c(
  "LT1"             = "#FF9300",
  "LT2"             = "#0000FF",
  "A26f"            = "#996633",
  "A18j"            = "#996633",
  "A01c"            = "#996633",
  "A14a"            = "#996633",
  "A09l (DnB)"      = "#FF2600",
  "A27k"            = "#FF66CC",
  "A02o (Wave-1)"   = "#800000",
  "TePn05"          = "#400040",
  "Goro"            = "#2C3E50",
  "A23g"            = "#FFD479",
  "A05q"            = "#FFD479",
  "A09c (Basin-4)"  = "#008F00",
  "A09g (Basin-3)"  = "#008F00",
  "A09a (Basin-2)"  = "#008F00",
  "chordotonal"     = "#004080",
  "chordotonals"    = "#004080",
  "mdIV"            = "#004080"
)

get_colors <- function(labels) {
  cols <- neuron_color_map[labels]
  if (any(is.na(cols))) {
    message("No color defined for: ",
            paste(unique(labels[is.na(cols)]), collapse = ", "),
            " -- defaulting to black.")
    cols[is.na(cols)] <- "black"
  }
  unname(cols)
}

# ---- 1. Load & clean --------------------------------------------------

raw <- read_excel("6-SRC_Proportion_downstream_summary_woutA19fl.xlsx",
                   sheet = "Proportion_downstream")

# Column headers have stray leading/trailing spaces (" LT1 ", " A26f ", ...)
names(raw) <- str_trim(names(raw))

# Source_neuron_name also has stray spaces in some rows (" A05q " etc.)
raw <- raw %>%
  mutate(Source_neuron_name = str_trim(Source_neuron_name))

target_cols <- setdiff(names(raw),
                        c("Source_neuron_CATMAID", "Source_neuron_name",
                          "Source_neuron_Segment", "Source_neuron_Side"))

# Some target columns can come in as character instead of numeric -- e.g.
# a literal "NA" text string (rather than a true blank/NA) anywhere in a
# column is enough to make readxl treat that whole column as character,
# which then fails when pivot_longer() tries to combine it with the
# numeric columns. Report exactly which cells this affects, then coerce.
for (col in target_cols) {
  if (!is.numeric(raw[[col]])) {
    bad_vals <- raw[[col]][is.na(suppressWarnings(as.numeric(raw[[col]]))) &
                              !is.na(raw[[col]])]
    if (length(bad_vals) > 0) {
      message("Column '", col, "' has non-numeric entries: ",
              paste(unique(bad_vals), collapse = ", "),
              " (", length(bad_vals), " cell(s)) -- coercing to NA.")
    }
  }
}
raw <- raw %>%
  mutate(across(all_of(target_cols), ~ suppressWarnings(as.numeric(.x))))

# Treat any missing proportion as 0 (no recorded connection) -- done here,
# numerically, on the whole pipeline, rather than relying on a manual edit
# in the Excel file (typing "0" into a cell by hand risks turning the
# column into character again if it's entered as text instead of a number,
# which silently breaks every ggplot/pheatmap fill downstream).
n_na <- sum(is.na(raw[target_cols]))
if (n_na > 0) {
  message(n_na, " missing target-column cell(s) are being treated as 0 ",
          "(no recorded connection).")
}
raw <- raw %>%
  mutate(across(all_of(target_cols), ~ replace_na(.x, 0)))

# Anatomical ordering: T1,T2,T3 before A1..A9; L before R
segment_levels <- c(paste0("T", 1:3), paste0("A", 1:9))
raw <- raw %>%
  mutate(
    Source_neuron_Segment = factor(Source_neuron_Segment, levels = segment_levels),
    Source_neuron_Side    = factor(Source_neuron_Side, levels = c("L", "R")),
    seg_side_label = paste0(Source_neuron_Segment, Source_neuron_Side)
  ) %>%
  arrange(Source_neuron_name, Source_neuron_Segment, Source_neuron_Side)

# Long format -- the workhorse for both ggplot2 and summarizing
long <- raw %>%
  pivot_longer(all_of(target_cols), names_to = "Target", values_to = "Proportion")

# Arcsine-square-root transform: standard variance-stabilizing transform
# for proportion/percentage data, useful here because the raw proportions
# span a wide range (near-0 up to ~20%+) and small differences among the
# low values would otherwise be visually swamped by the few large ones.
# Proportion is a percentage (0-100), so convert to a 0-1 fraction first.
long <- long %>%
  mutate(Proportion_arcsin = asin(sqrt(Proportion / 100)))



# ---- 1. Group-level summary matrix (median proportion per group pair) -
# Collapses the 148 individual rows down to one row per source-neuron
# GROUP (e.g. all mdIV_a*l/r) -- the natural "who talks to whom" summary.
# Uses the median (not mean) across segmental instances, matching
# connectivity_fan_groups_v10.R's summarise_group().
#
# IMPORTANT: the median is taken on the RAW percentage first, and the
# arcsin(sqrt(.)) transform applied AFTER, to that one summary value --
# not the other way around. Median and a nonlinear transform only commute
# when a group has an odd number of values (the median is then a single
# actual observation, and transforming a single point trivially
# commutes). For an EVEN-sized group, the conventional median is the
# average of the two middle values, and that averaging step does not
# commute with a nonlinear function: arcsin(sqrt((a+b)/2)) is not equal
# to (arcsin(sqrt(a)) + arcsin(sqrt(b))) / 2. Many source/target groups
# here have an even number of segmental instances, so taking the median
# after transforming (as an earlier version of this script did) would
# give a different, incorrect answer.

group_matrix <- long %>%
  group_by(Source_neuron_name, Target) %>%
  summarise(median_prop_raw = median(Proportion, na.rm = TRUE), .groups = "drop") %>%
  mutate(median_prop = asin(sqrt(median_prop_raw / 100))) %>%
  select(-median_prop_raw) %>%
  pivot_wider(names_from = Target, values_from = median_prop) %>%
  column_to_rownames("Source_neuron_name") %>%
  as.matrix()

# 1a. Clustered heatmap -- reveals which source groups share a similar
#     downstream "fingerprint", and which target groups tend to be hit
#     together.
#
# The dendrograms are still built from the data (same distance/linkage
# pheatmap would use internally), but branches are then ROTATED (not
# restructured) so that, wherever the clustering allows it, groups fall
# into the circuit-hierarchy order: sensory (mdIV, chordotonal) -> basins
# -> other interneurons -> preLTs -> motoneurons (LT1/LT2), i.e. low rank
# = top/left, high rank = bottom/right. This only swaps sibling branches
# at existing merge points -- it can't force perfect separation between
# groups if the data itself clusters them differently, but it will get as
# close as the tree topology allows.


hierarchy_rank <- c(
  "mdIV" = 1, "chordotonal" = 1, "chordotonals" = 1,
  "A09a (Basin-2)" = 2, "A09a (Basin-2)" = 2, "A09g (Basin-3)" = 2, "A09c (Basin-4)" = 2,
  "A05q" = 3, "A23g" = 3, "Goro" = 4, "TePn05" = 4, "A02o (Wave-1)" = 4, 
  "A14a" = 6, "A01c" = 6, "A18j" = 6, "A26f" = 6, "A09l (DnB)" = 4, "A27k" = 5,
  "LT1" = 7, "LT2" = 7
)
get_rank <- function(labels) {
  r <- hierarchy_rank[labels]
  if (any(is.na(r))) {
    message("No hierarchy rank defined for: ",
            paste(unique(labels[is.na(r)]), collapse = ", "),
            " -- leaving their position to the clustering (mid-rank fallback).")
    r[is.na(r)] <- 3
  }
  r
}

reorder_by_hierarchy <- function(mat, margin = c("row", "col")) {
  margin <- match.arg(margin)
  m <- if (margin == "row") mat else t(mat)
  hc <- hclust(dist(m, method = "euclidean"), method = "average")
  wts <- get_rank(hc$labels)
  as.hclust(reorder(as.dendrogram(hc), wts, agglo.FUN = mean))
}

row_hc <- reorder_by_hierarchy(group_matrix, "row")
col_hc <- reorder_by_hierarchy(group_matrix, "col")

# Build pheatmap silently (no drawing yet) so we can recolor the row/column
# name labels before writing the final PDF. pheatmap always lays out labels
# in the post-clustering (dendrogram) order, given by tree_row$order /
# tree_col$order, so the color vectors below must follow that same order,
# not the original row/column order of group_matrix.
p_heat <- pheatmap(group_matrix,
                    color = colorRampPalette(c("white", "orange", "firebrick"))(100),
                    cluster_rows = row_hc,
                    cluster_cols = col_hc,
                    main = "Median arcsin(sqrt(proportion)), by source group\n(rotated toward sensory -> motor hierarchy)",
                    silent = TRUE)

row_order  <- p_heat$tree_row$order
col_order  <- p_heat$tree_col$order
row_colors <- get_colors(rownames(group_matrix)[row_order])
col_colors <- get_colors(colnames(group_matrix)[col_order])

# Recolor the row_names / col_names text grobs in place. This relies on
# pheatmap's internal gtable layout names ("row_names", "col_names"), which
# have been stable across pheatmap versions; if a future pheatmap version
# renames them, this step is skipped with a warning rather than erroring
# out, and the heatmap still saves with default black labels.
recolor_pheatmap_labels <- function(pheatmap_obj, row_colors, col_colors) {
  g <- pheatmap_obj$gtable
  row_idx <- which(g$layout$name == "row_names")
  col_idx <- which(g$layout$name == "col_names")
  if (length(row_idx) == 1) {
    g$grobs[[row_idx]]$gp$col <- row_colors
  } else {
    warning("Could not find 'row_names' grob -- row labels left black.")
  }
  if (length(col_idx) == 1) {
    g$grobs[[col_idx]]$gp$col <- col_colors
  } else {
    warning("Could not find 'col_names' grob -- column labels left black.")
  }
  g
}

g_colored <- recolor_pheatmap_labels(p_heat, row_colors, col_colors)

pdf("group_level_heatmap_outA19fl.pdf", width = 10, height = 8)
grid.newpage()
grid.draw(g_colored)
dev.off()

# 1b. Chord diagram -- compact circular summary of the same group-level
#     matrix; edge thickness = median transformed proportion. Sector colors
#     (and, by extension, the connections leaving each sector, which take
#     their source sector's color by default) follow the same neuron color
#     scheme used above.
# NOTE: values are on the arcsin(sqrt(p)) scale (range ~0-1.57 radians,
# not 0-100%), so the filter threshold below is on that scale -- tune it
# to taste depending on how cluttered/sparse the diagram looks.
chord_df <- group_matrix %>%
  as.data.frame() %>%
  rownames_to_column("source") %>%
  pivot_longer(-source, names_to = "target", values_to = "value") %>%
  filter(value > 0.05)

chord_sectors <- union(chord_df$source, chord_df$target)
chord_grid_col <- setNames(get_colors(chord_sectors), chord_sectors)

pdf("chord_diagram.pdf", width = 9, height = 9)
circos.clear()
chordDiagram(chord_df, grid.col = chord_grid_col, transparency = 0.3,
             annotationTrack = c("grid"),
             preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector <- CELL_META$sector.index
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector,
               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.7,
               col = chord_grid_col[sector])
}, bg.border = NA)
dev.off()
circos.clear()

