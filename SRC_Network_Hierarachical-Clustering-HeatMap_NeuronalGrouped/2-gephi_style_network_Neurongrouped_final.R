## ---------------------------------------------------------------------
## Gephi-style conceptual network plot of the group-level connectivity
## matrix (same group_matrix / chord_df data as section 1a/1b of
## 1-analyze_connectivity_matrix_final.R), styled after force-directed
## network figures
## (colored node clusters, node size ~ connectivity strength, curved
## directed edges colored by their source node).
##
## This script is self-contained: it repeats the minimal data prep needed
## to rebuild group_matrix, then builds the network. If you've already
## run 1-analyze_connectivity_matrix_final.R in the same R session,
## you can skip straight to "PART B" below and reuse its group_matrix /
## neuron_color_map directly.
## ---------------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(shadowtext)

# =========================================================================
# PART A -- rebuild group_matrix (identical logic to the v7 script)
# =========================================================================

raw <- read_excel("6-SRC_Proportion_downstream_summary_woutA19fl.xlsx",
                  sheet = "Proportion_downstream")
names(raw) <- str_trim(names(raw))
raw <- raw %>% mutate(Source_neuron_name = str_trim(Source_neuron_name))

target_cols <- setdiff(names(raw),
                       c("Source_neuron_CATMAID", "Source_neuron_name",
                         "Source_neuron_Segment", "Source_neuron_Side"))

raw <- raw %>%
  mutate(across(all_of(target_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(target_cols), ~ replace_na(.x, 0)))

long <- raw %>%
  pivot_longer(all_of(target_cols), names_to = "Target", values_to = "Proportion") %>%
  mutate(Proportion_arcsin = asin(sqrt(Proportion / 100)))

# Median is taken on the RAW percentage first, then arcsin(sqrt(.))
# applied after -- see analyze_connectivity_matrix_v7_colored.R for why
# (median and a nonlinear transform don't commute for even-sized groups,
# which many source/target groups here are).
group_matrix <- long %>%
  group_by(Source_neuron_name, Target) %>%
  summarise(median_prop_raw = median(Proportion, na.rm = TRUE), .groups = "drop") %>%
  mutate(median_prop = asin(sqrt(median_prop_raw / 100))) %>%
  select(-median_prop_raw) %>%
  pivot_wider(names_from = Target, values_from = median_prop) %>%
  column_to_rownames("Source_neuron_name") %>%
  as.matrix()

# Same fixed color scheme as v7 (both "chordotonal"/"chordotonals" keys kept
# so the color applies whether the name shows up as a source row or a
# target column).
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
            " -- defaulting to grey60.")
    cols[is.na(cols)] <- "grey60"
  }
  unname(cols)
}

# =========================================================================
# PART B -- build the network graph
# =========================================================================

# Edge-weight threshold on the arcsin(sqrt(p)) scale, same scale as the
# chord diagram in v7 -- tune this if the plot looks too dense/sparse.
EDGE_THRESHOLD <- 0.000

edges <- group_matrix %>%
  as.data.frame() %>%
  rownames_to_column("source") %>%
  pivot_longer(-source, names_to = "target", values_to = "weight") %>%
  filter(weight > EDGE_THRESHOLD, source != target)

# One node per neuron/group appearing on EITHER side (source or target) --
# this is a general directed graph, not a strict bipartite one, since
# several groups (e.g. A09l (DnB)) act as both a source and a downstream
# target elsewhere in the matrix.
all_names <- union(edges$source, edges$target)

# Hierarchy tier -- same 5-tier sensory -> motor scheme used to rotate the
# v7 heatmap dendrogram, reused here to pin vertical position: tier 1
# (sensory: mdIV, chordotonal) is forced to the top of the plot, tier 5
# (LT1/LT2) to the bottom, with everything else free to fall in between
# tiers 2-4. Within a tier, horizontal (x) position is still determined by
# force-directed layout, so connected/similar nodes still cluster together
# left-right -- only the vertical axis is constrained.
hierarchy_tier <- c(
  "mdIV" = 1, "chordotonal" = 1, "chordotonals" = 1,
  "A09a (Basin-2)" = 2, "A09a_(Basin-2)" = 2, "A09g (Basin-3)" = 2, "A09c (Basin-4)" = 2,
  "A05q" = 3, "A23g" = 3, "Goro" = 4, "TePn05" = 4, "A02o (Wave-1)" = 4,
   "A14a" = 6, "A01c" = 6, "A18j" = 6, "A26f" = 6, "A09l (DnB)" = 4, "A27k" = 5,
  "LT1" = 7, "LT2" = 7
)

get_tier <- function(labels) {
  t <- hierarchy_tier[labels]
  if (any(is.na(t))) {
    message("No hierarchy tier defined for: ",
            paste(unique(labels[is.na(t)]), collapse = ", "),
            " -- defaulting to mid-tier (3).")
    t[is.na(t)] <- 3
  }
  unname(t)
}

# Neuron type -> node shape: triangle for sensory, hexagon (flat edge on
# top) for interneurons, circle for motoneurons.
neuron_type_map <- c(
  "mdIV" = "sensory", "chordotonal" = "sensory", "chordotonals" = "sensory",
  "A09a (Basin-2)" = "interneuron", "A09a_(Basin-2)" = "interneuron",
  "A09g (Basin-3)" = "interneuron", "A09c (Basin-4)" = "interneuron",
  "A05q" = "interneuron", "A23g" = "interneuron",
  "Goro" = "interneuron", "TePn05" = "interneuron",
  "A02o (Wave-1)" = "interneuron", "A09l (DnB)" = "interneuron",
  "A27k" = "interneuron", "A26f" = "interneuron",
  "A18j" = "interneuron", "A01c" = "interneuron", "A14a" = "interneuron",
  "LT1" = "motoneuron", "LT2" = "motoneuron"
)

get_type <- function(labels) {
  ty <- neuron_type_map[labels]
  if (any(is.na(ty))) {
    message("No neuron type defined for: ",
            paste(unique(labels[is.na(ty)]), collapse = ", "),
            " -- defaulting to interneuron (hexagon).")
    ty[is.na(ty)] <- "interneuron"
  }
  unname(ty)
}

nodes <- tibble(name = all_names) %>%
  mutate(color = get_colors(name),
         tier  = get_tier(name),
         type  = get_type(name))

g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

set.seed(42)  # force-directed layout is stochastic -- fix for reproducibility

# Build a base force-directed layout first (for natural left-right
# clustering of connected nodes), then override the y-coordinate with the
# fixed tier, and re-space x within each tier so nodes never overlap --
# relative left-right order within a tier is preserved from the
# force-directed pass, only exact spacing is evened out.
lay <- create_layout(g, layout = "fr")

lay <- lay %>%
  group_by(tier) %>%
  mutate(
    x = {
      n_in_tier <- n()
      spacing <- 1.6
      positions <- seq(-(n_in_tier - 1) / 2, (n_in_tier - 1) / 2, length.out = n_in_tier) * spacing
      positions[rank(x, ties.method = "first")]
    }
  ) %>%
  ungroup() %>%
  mutate(y = -tier)  # tier 1 (sensory) -> least negative -> plotted at top;
# tier 5 (LT1/LT2) -> most negative -> plotted at bottom

# ---- Node shapes -------------------------------------------------------
# ggplot2's built-in point shapes (pch 0-25) don't include a flat-top
# hexagon, so all three node shapes are drawn as explicit polygons at each
# node's (x, y) instead of via geom_point()'s shape aesthetic -- this
# keeps triangle/hexagon/circle visually consistent (same drawing method,
# same outline style) and gives exact control over the hexagon's
# orientation (flat edge on top, as requested).
#
# NODE_RADIUS is in the same units as the layout's x/y coordinates (which
# span roughly -5..5 in x within a tier and -5..-1 in y across tiers) --
# tune this if nodes look too small/large or start overlapping once you
# see the actual plot.
NODE_RADIUS <- 0.28

make_polygon_vertices <- function(cx, cy, radius, shape) {
  angles_deg <- switch(shape,
                       "triangle" = c(90, 210, 330),           # point-up triangle
                       "hexagon"  = seq(0, 300, by = 60),       # flat edge on top (no vertex at 90 deg)
                       "circle"   = seq(0, 350, by = 9),        # 40-gon approximating a circle
                       stop("Unknown node shape: ", shape)
  )
  rad <- angles_deg * pi / 180
  tibble(px = cx + radius * cos(rad), py = cy + radius * sin(rad))
}

node_polygons <- lay %>%
  as_tibble() %>%
  select(name, x, y, color, type) %>%
  mutate(shape = case_when(
    type == "sensory"     ~ "triangle",
    type == "interneuron" ~ "hexagon",
    type == "motoneuron"  ~ "circle",
    TRUE                  ~ "circle"
  )) %>%
  rowwise() %>%
  mutate(poly = list(make_polygon_vertices(x, y, NODE_RADIUS, shape))) %>%
  ungroup() %>%
  unnest(poly)

# ---- Edge-width legend (discrete intervals, not node-size) ------------
# Connectivity strength is now shown purely by edge (line) width, binned
# into intervals that span the full range of the plotted edge weights, so
# the legend reads as a small set of labeled line-width swatches rather
# than a continuous node-size gradient.
N_STRENGTH_BINS <- 5
weight_breaks <- pretty(range(edges$weight), n = N_STRENGTH_BINS)
weight_breaks[1] <- min(weight_breaks[1], min(edges$weight))
weight_breaks[length(weight_breaks)] <- max(weight_breaks[length(weight_breaks)], max(edges$weight))

edges <- edges %>%
  mutate(weight_bin = cut(weight, breaks = weight_breaks, include.lowest = TRUE, dig.lab = 3))

width_values <- setNames(seq(0.4, 3.5, length.out = nlevels(edges$weight_bin)),
                         levels(edges$weight_bin))

g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)  # rebuild with weight_bin column

p_gephi <- ggraph(g, layout = "manual", x = lay$x, y = lay$y) +
  geom_edge_fan(
    aes(width = weight_bin, alpha = weight, color = node1.color),
    arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
    end_cap = circle(NODE_RADIUS * 10, "mm"),
    start_cap = circle(NODE_RADIUS * 10, "mm"),
    strength = 0.6
  ) +
  scale_edge_width_manual(values = width_values, name = "Connectivity\nstrength", drop = FALSE) +
  scale_edge_alpha(range = c(0.3, 0.85), guide = "none") +
  scale_edge_color_identity() +
  geom_polygon(data = node_polygons, aes(x = px, y = py, group = name, fill = color),
               color = "black", linewidth = 0.4) +
  scale_fill_identity() +
  geom_shadowtext(data = as_tibble(lay) %>% select(name, x, y),
                  aes(x = x, y = y, label = name),
                  color = "black", bg.colour = "white", bg.r = 0.12,
                  size = 5, fontface = "bold", lineheight = 0.85) +
  theme_void() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold")) +
  guides(edge_width = guide_legend(override.aes = list(alpha = 1))) +
  labs(title = "Conceptual connectivity network (group-level, arcsin(sqrt(p)) weighted)",
       subtitle = paste0("Edges above threshold ", EDGE_THRESHOLD,
                         "; sensory neurons (mdIV, chordotonal) pinned top, LT1/LT2 pinned bottom\n",
                         "Shapes: triangle = sensory, hexagon = interneuron, circle = motoneuron; ",
                         "colors = neuron group; line width = connectivity strength")) +
  coord_fixed()

ggsave("gephi_style_network_2_woutA19fl.pdf", p_gephi, width = 11, height = 9)
ggsave("gephi_style_network_2_woutA19fl.png", p_gephi, width = 11, height = 9, dpi = 300)
# Requires the 'svglite' package (install.packages("svglite")) -- ggsave()
# uses it automatically for .svg output when available.
ggsave("gephi_style_network_2_woutA19fl.svg", p_gephi, width = 11, height = 9)

