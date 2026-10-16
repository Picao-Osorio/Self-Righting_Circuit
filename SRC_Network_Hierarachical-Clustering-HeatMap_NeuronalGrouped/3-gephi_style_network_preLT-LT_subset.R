## ---------------------------------------------------------------------
## Gephi-style conceptual network plot -- LT-PATHWAY SUBSET ONLY
## Same data/pipeline as gephi_style_network.R, restricted to:
##   LT1, LT2, A26f, A18j, A01c, A14a, A09l (DnB), A27k
## (i.e. the preLT + prepreLT groups from connectivity_fan_groups_v10.R,
## combined into one network instead of separate per-group fan diagrams.)
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
# PART A -- rebuild group_matrix (identical logic to the full script)
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

# Median taken on the RAW percentage first, then arcsin(sqrt(.)) applied
# after (median and a nonlinear transform don't commute for even-sized
# groups -- see analyze_connectivity_matrix_v7_colored.R for the full
# explanation).
group_matrix <- long %>%
  group_by(Source_neuron_name, Target) %>%
  summarise(median_prop_raw = median(Proportion, na.rm = TRUE), .groups = "drop") %>%
  mutate(median_prop = asin(sqrt(median_prop_raw / 100))) %>%
  select(-median_prop_raw) %>%
  pivot_wider(names_from = Target, values_from = median_prop) %>%
  column_to_rownames("Source_neuron_name") %>%
  as.matrix()

# ---- Restrict to the requested subset ----------------------------------
# Keep only these neurons on BOTH axes (as source rows and as target
# columns) -- anything connecting only to a neuron outside this list is
# dropped along with it.
subset_names <- c("LT1", "LT2", "A26f", "A18j", "A01c", "A14a", "A09l (DnB)", "A27k")

missing_from_matrix <- setdiff(subset_names, union(rownames(group_matrix), colnames(group_matrix)))
if (length(missing_from_matrix) > 0) {
  warning("Not found in group_matrix (check spelling against the Excel headers): ",
          paste(missing_from_matrix, collapse = ", "))
}

group_matrix <- group_matrix[rownames(group_matrix) %in% subset_names,
                              colnames(group_matrix) %in% subset_names,
                              drop = FALSE]

# Same fixed color scheme as the full script (only the entries relevant to
# this subset are used, the rest are simply unused).
neuron_color_map <- c(
  "LT1"             = "#FF9300",
  "LT2"             = "#0000FF",
  "A26f"            = "#996633",
  "A18j"            = "#996633",
  "A01c"            = "#996633",
  "A14a"            = "#996633",
  "A09l (DnB)"      = "#FF2600",
  "A27k"            = "#FF66CC",
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

EDGE_THRESHOLD <- 0.001

edges <- group_matrix %>%
  as.data.frame() %>%
  rownames_to_column("source") %>%
  pivot_longer(-source, names_to = "target", values_to = "weight") %>%
  filter(weight > EDGE_THRESHOLD, source != target)

all_names <- union(edges$source, edges$target)

# Hierarchy tier -- only the tiers relevant to this subset : A27k/A09l (DnB) sit at tier 3-4, the five preLTs
# at tier 4, LT1/LT2 pinned to tier 5 (bottom).

hierarchy_tier <- c(
   "A09l (DnB)" = 1, "A27k" = 1,
   "A26f" = 2, "A18j" = 2, "A01c" = 2, "A14a" = 2, 
  "LT1" = 3, "LT2" = 3
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

# Neuron type -> node shape: this subset has no sensory neurons, so only
# hexagon (interneuron) and circle (motoneuron) shapes appear.
neuron_type_map <- c(
  "A27k" = "interneuron",  "A26f" = "interneuron",
  "A18j" = "interneuron", "A01c" = "interneuron", "A14a" = "interneuron",
   "A09l (DnB)" = "interneuron",
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

set.seed(42)

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
  mutate(y = -tier)

NODE_RADIUS <- 0.28

make_polygon_vertices <- function(cx, cy, radius, shape) {
  angles_deg <- switch(shape,
    "triangle" = c(90, 210, 330),
    "hexagon"  = seq(0, 300, by = 60),
    "circle"   = seq(0, 350, by = 9),
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

N_STRENGTH_BINS <- 5
weight_breaks <- pretty(range(edges$weight), n = N_STRENGTH_BINS)
weight_breaks[1] <- min(weight_breaks[1], min(edges$weight))
weight_breaks[length(weight_breaks)] <- max(weight_breaks[length(weight_breaks)], max(edges$weight))

edges <- edges %>%
  mutate(weight_bin = cut(weight, breaks = weight_breaks, include.lowest = TRUE, dig.lab = 3))

width_values <- setNames(seq(0.4, 3.5, length.out = nlevels(edges$weight_bin)),
                          levels(edges$weight_bin))

g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

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
  labs(title = "Conceptual connectivity network -- LT pathway subset",
       subtitle = paste0("Edges above threshold ", EDGE_THRESHOLD,
                          "; LT1/LT2 pinned to bottom\n",
                          "Shapes: hexagon = interneuron, circle = motoneuron; ",
                          "colors = neuron group; line width = connectivity strength")) +
  coord_fixed()

ggsave("gephi_style_network_LT_subset_2.pdf", p_gephi, width = 8, height = 7)
ggsave("gephi_style_network_LT_subset_2.png", p_gephi, width = 8, height = 7, dpi = 300)
# Requires the 'svglite' package (install.packages("svglite"))
ggsave("gephi_style_network_LT_subset.svg", p_gephi, width = 8, height = 7)
