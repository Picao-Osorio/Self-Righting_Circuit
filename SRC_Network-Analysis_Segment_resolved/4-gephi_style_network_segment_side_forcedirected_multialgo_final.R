## ---------------------------------------------------------------------
## Fully segment/side-resolved Gephi-style conceptual network plots --
## FRUCHTERMAN-REINGOLD force-directed layout, circuit MINUS A19f/A19l only.
##

##
## Layout: every plot is driven SOLELY by the graph's own connectivity, via
## Fruchterman-Reingold (ggraph layout = "fr" -> igraph::layout_with_fr(),
## which is weight-aware but ignores edge direction). No hierarchy is imposed.
## The only post-processing is a pure rotation/reflection (orientation), which
## never changes any node's position relative to any other.
##
## Dimensions:
##   - resolution: side x side (fully resolved) AND seg x seg (L/R averaged)
##   - weight: raw #Synapses AND arcsin(sqrt(proportion))
##   - circuit: MINUS A19f/A19l (nodes whose Group is A19f or A19l, and every
##     edge touching them, are dropped before layout)
## => 2 resolutions x 2 weight kinds = 4 plots (each saved as .pdf/.png/.svg).
##
## SOURCE NEURONS ONLY: only nodes that are SOURCE
## neurons/instances are laid out and plotted. Source neurons are those of the
## workbook's Filtered_downstream sheet (source_neuron_* columns); checked
## against the data prep script's output, they are exactly the instances that
## have a ROW in the connectivity matrix. Instances that only ever appear as a
## downstream TARGET (a matrix column with no matching row -- e.g. A09l (DnB)
## at T3 and A8) are dropped, together with every edge ending on them, BEFORE
## the layout is computed, so they do not influence the placement of the
## remaining nodes. EXCEPTION: the LT1 and LT2 motoneurons (A1-A7) are
## targets only but are always kept (SOURCE_EXCEPTION_GROUPS). Set
## SOURCE_NODES_ONLY <- FALSE to plot all nodes. A source whose only
## above-threshold edges went to dropped nodes has no edges left and
## therefore does not appear either.
## ---------------------------------------------------------------------

required_pkgs <- c("igraph", "ggraph", "tidygraph", "ggplot2", "shadowtext")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
invisible(lapply(required_pkgs, library, character.only = TRUE))

source("0-segment_side_data_prep_final.R")

# ---- Build a unified node-metadata table from a row_meta/col_meta table -
# (Row_label_side/Row_label_seg/Col_label_side/Col_label_seg, Group, Seg,
# [Side]) -- adds Side = "none" when the table has no Side column (the
# seg-only metadata tables).

meta_to_nodes <- function(meta, name_col) {
  meta <- meta %>% rename(name = !!name_col)
  if (!"Side" %in% names(meta)) meta$Side <- "none"
  meta %>% select(name, Group, Seg, Side)
}

# Display-only label simplification (unchanged from the previous script).
display_name_map <- c(
  "A09l (DnB)"     = "DnB",
  "A02o (Wave-1)"  = "Wave-1",
  "A09c (Basin-4)" = "Basin-4",
  "A09g (Basin-3)" = "Basin-3",
  "A09a (Basin-2)" = "Basin-2",
  "chordotonals"   = "Ch"
)

get_display_group <- function(groups) {
  simplified <- unname(display_name_map[groups])
  ifelse(is.na(simplified), groups, simplified)
}

# Groups to drop for the "minus A19" circuit version -- 
# request, same pair excluded throughout the rest of this project
# (structural_network_helpers.R's EXCLUDE_A19).
EXCLUDE_A19 <- c("A19f", "A19l")

SOURCE_EXCEPTION_GROUPS <- c("LT1", "LT2")  # motoneurons kept although they are targets only
SOURCE_NODES_ONLY <- TRUE  # TRUE = plot only source neurons (nodes with a row in the matrix); FALSE = old behaviour

ALGO_LABEL <- c(fr = "Fruchterman-Reingold")

# ---- Shared network-building + plotting routine ------------------------

build_gephi_network <- function(mat, row_meta, row_name_col, col_meta, col_name_col,
                                 weight_kind = c("raw", "arcsin"),
                                 edge_threshold, out_prefix,
                                 layout_algo = "fr",
                                 exclude_groups = NULL,
                                 fig_width, fig_height) {
  weight_kind <- match.arg(weight_kind)

  # check.names = FALSE is essential here: several canonical Group names
  # contain parentheses/spaces (e.g. "A09a (Basin-2)"), and as.data.frame()'s
  # default check.names = TRUE runs make.names() on the matrix's column
  # names, silently mangling them. Row names are untouched by
  # rownames_to_column(), so only TARGET nodes would be affected.
  edges <- mat %>%
    as.data.frame(check.names = FALSE) %>%
    rownames_to_column("source") %>%
    pivot_longer(-source, names_to = "target", values_to = "weight") %>%
    filter(!is.na(weight), weight > edge_threshold)

  if (nrow(edges) == 0) {
    warning("No edges above threshold ", edge_threshold, " for ", out_prefix, " -- skipping.")
    return(invisible(NULL))
  }

  # ---- Source neurons only: drop target-only instances (no row in the
  # matrix) and every edge ending on them, before nodes/layout are built. ----
  if (SOURCE_NODES_ONLY) {
    exception_names <- col_meta[[col_name_col]][col_meta$Group %in% SOURCE_EXCEPTION_GROUPS]
    source_names <- union(unique(row_meta[[row_name_col]]), exception_names)
    n_target_only <- length(setdiff(unique(edges$target), source_names))
    edges <- edges %>% filter(source %in% source_names, target %in% source_names)
    message("    source-only filter: dropped ", n_target_only, " target-only node(s) (kept: ", paste(SOURCE_EXCEPTION_GROUPS, collapse = "/"), ") for ", out_prefix, ".")
    if (nrow(edges) == 0) {
      warning("No edges left among source neurons for ", out_prefix, " -- skipping.")
      return(invisible(NULL))
    }
  }

  all_nodes_meta <- bind_rows(
    meta_to_nodes(row_meta, row_name_col),
    meta_to_nodes(col_meta, col_name_col)
  ) %>% distinct(name, .keep_all = TRUE)

  node_names <- union(unique(edges$source), unique(edges$target))
  nodes <- tibble(name = node_names) %>%
    left_join(all_nodes_meta, by = "name") %>%
    mutate(color = get_colors_by_group(Group),
           type  = get_type(Group))

  if (any(is.na(nodes$Group))) {
    stop("Internal error: node(s) with no metadata match: ", paste(nodes$name[is.na(nodes$Group)], collapse = ", "))
  }

  # ---- Optional circuit exclusion (A19f/A19l for the "minus A19"
  # version): drop AFTER building the node table so the same
  # node-identification logic is reused for both circuit versions; this
  # also drops every edge touching an excluded node, rather than leaving
  # zero-degree remnants behind (same approach as
  # structural_network_helpers.R's build_structural_graph()). ----
  if (!is.null(exclude_groups) && length(exclude_groups) > 0) {
    keep_names <- nodes$name[!(nodes$Group %in% exclude_groups)]
    nodes <- nodes %>% filter(name %in% keep_names)
    edges <- edges %>% filter(source %in% keep_names, target %in% keep_names)
    if (nrow(edges) == 0) {
      warning("No edges remain for ", out_prefix, " after excluding ",
              paste(exclude_groups, collapse = "/"), " -- skipping.")
      return(invisible(NULL))
    }
  }

  nodes$segi <- seg_index(nodes$Seg)

  # ---- Display label (plot text only -- see display_name_map above). ----
  nodes$display_label <- paste0(
    get_display_group(nodes$Group), "_",
    ifelse(nodes$Side == "none", nodes$Seg, paste0(nodes$Seg, nodes$Side))
  )

  g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

  set.seed(42)
  lay <- create_layout(g, layout = "fr")

  # Orientation. Fruchterman-Reingold has no notion of "up" or "left", so
  # its output orientation is arbitrary. This block only ROTATES/REFLECTS the
  # coordinates it already produced -- every pairwise distance and angle is
  # preserved, so the force algorithm alone still determines the topology;
  # it just resolves orientation for visual consistency (like orienting a
  # map to north). The whole layout is rotated about its centroid so the
  # vector from the sensory centroid to the motoneuron centroid points
  # straight down.
  sensory_idx <- lay$type == "sensory"
  motor_idx   <- lay$type == "motoneuron"
  if (any(sensory_idx) && any(motor_idx)) {
    sx <- mean(lay$x[sensory_idx]); sy <- mean(lay$y[sensory_idx])
    mx <- mean(lay$x[motor_idx]);   my <- mean(lay$y[motor_idx])
    vx <- mx - sx; vy <- my - sy
    if (vx != 0 || vy != 0) {
      theta <- atan2(vy, vx)
      alpha <- -pi / 2 - theta  # rotate so the sensory->motor vector points straight down (negative y)
      cx <- mean(lay$x); cy <- mean(lay$y)
      xs <- lay$x - cx; ys <- lay$y - cy
      lay$x <- xs * cos(alpha) - ys * sin(alpha) + cx
      lay$y <- xs * sin(alpha) + ys * cos(alpha) + cy
    }
  }

  # Left-right (axial segment, T1->T3/A1->A9 increasing left to right): a
  # mirror (x-sign flip only, y untouched) rather than a further rotation,
  # since the up-down orientation fixed just above must be preserved --
  # mirroring x alone leaves y, and therefore the sensory/motor
  # orientation, exactly as set.
  if (length(unique(lay$segi)) > 1) {
    x_seg_cor <- suppressWarnings(stats::cor(lay$x, lay$segi))
    if (!is.na(x_seg_cor) && x_seg_cor < 0) {
      lay$x <- -lay$x
    }
  }

  # ---- Node shapes (triangle=sensory, hexagon=interneuron, circle=motoneuron) ----
  # NODE_RADIUS is derived as a fraction of the layout's own coordinate
  # spread (FR/KK coordinates have no fixed scale) to keep node markers
  # sensibly sized across circuit versions / algorithms of different sizes.
  NODE_RADIUS <- max(0.15, diff(range(c(lay$x, lay$y))) * 0.012)
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

  # ---- Edge-width legend: discrete intervals spanning the observed range ----
  # Raw #Synapses uses a fixed, hand-picked set of breaks (WEIGHT_BREAKS_RAW);
  # arcsin(sqrt(proportion)) uses automatic, evenly-spaced pretty() breaks.
  if (weight_kind == "raw") {
    weight_breaks <- WEIGHT_BREAKS_RAW
  } else {
    N_STRENGTH_BINS <- 5
    weight_breaks <- pretty(range(edges$weight), n = N_STRENGTH_BINS)
  }
  weight_breaks[1] <- min(weight_breaks[1], min(edges$weight))
  weight_breaks[length(weight_breaks)] <- max(weight_breaks[length(weight_breaks)], max(edges$weight))

  edges <- edges %>%
    mutate(weight_bin = cut(weight, breaks = weight_breaks, include.lowest = TRUE, dig.lab = 3))

  width_values <- setNames(seq(0.3, 3, length.out = nlevels(edges$weight_bin)),
                            levels(edges$weight_bin))

  # Edge opacity, binned the same way as width rather than mapped
  # continuously from raw weight (raw #Synapses is heavily right-skewed,
  # which compresses a continuous alpha scale toward its low end and makes
  # most arrows hard to see).
  alpha_values <- setNames(seq(0.5, 0.95, length.out = nlevels(edges$weight_bin)),
                            levels(edges$weight_bin))

  g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

  legend_name <- if (weight_kind == "raw") "#Synapses" else "arcsin(sqrt(p))"
  weight_desc <- if (weight_kind == "raw") "raw #Synapses" else "arcsin(sqrt(#Synapses / source's Total #Synapses))"

  circuit_desc <- if (!is.null(exclude_groups) && length(exclude_groups) > 0) {
    paste0("Circuit minus ", paste(exclude_groups, collapse = "/"))
  } else {
    "Full circuit"
  }

  layout_desc <- paste0(
    "No hierarchy imposed -- plain force-directed (", ALGO_LABEL[[layout_algo]], ") layout; ",
    "both position axes are determined solely by the graph's own connectivity. ", circuit_desc,
    if (SOURCE_NODES_ONLY) "; source neurons (+ LT1/LT2 motoneurons) only." else ".", ""
  )

  p_gephi <- ggraph(g, layout = "manual", x = lay$x, y = lay$y) +
    geom_edge_fan(
      aes(width = weight_bin, alpha = weight_bin, color = node1.color),
      arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
      end_cap = circle(NODE_RADIUS * 10, "mm"),
      start_cap = circle(NODE_RADIUS * 10, "mm"),
      strength = 0.5
    ) +
    scale_edge_width_manual(values = width_values, name = legend_name, drop = FALSE) +
    scale_edge_alpha_manual(values = alpha_values, guide = "none", drop = FALSE) +
    scale_edge_color_identity() +
    geom_polygon(data = node_polygons, aes(x = px, y = py, group = name, fill = color),
                 color = "black", linewidth = 0.3) +
    scale_fill_identity() +
    geom_shadowtext(data = as_tibble(lay) %>% select(name, display_label, x, y),
                    aes(x = x, y = y, label = display_label),
                    color = "black", bg.colour = "white", bg.r = 0.1,
                    size = 1.8, fontface = "bold", lineheight = 0.8) +
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold")) +
    guides(edge_width = guide_legend(override.aes = list(alpha = 1))) +
    labs(title = paste0("Conceptual connectivity network -- fully segment/side-resolved (", weight_desc,
                         "; ", ALGO_LABEL[[layout_algo]], " layout, no hierarchy)"),
         subtitle = paste0(layout_desc,
                            "\nShapes: triangle = sensory, hexagon = interneuron, circle = motoneuron; ",
                            "colors = neuron group; line width = ", legend_name)) +
    coord_fixed()

  ggsave(paste0(out_prefix, ".pdf"), p_gephi, width = fig_width, height = fig_height, limitsize = FALSE)
  ggsave(paste0(out_prefix, ".png"), p_gephi, width = fig_width, height = fig_height, dpi = 300, limitsize = FALSE)
  # Requires the 'svglite' package (install.packages("svglite"))
  ggsave(paste0(out_prefix, ".svg"), p_gephi, width = fig_width, height = fig_height, limitsize = FALSE)
  message("Saved ", out_prefix, ".{pdf,png,svg} (", nrow(nodes), " nodes, ", nrow(edges), " edges).")
}

# Edge-weight thresholds -- filter() keeps weight > threshold (strict), so
# EDGE_THRESHOLD_RAW <- 2 keeps only edges with MORE than 2 synapses
# (i.e. >=4). EDGE_THRESHOLD_ARCSIN is left at 0 (any nonzero connection
# shown), same as the original hybrid_6 script -- set it explicitly if an
# equivalent cutoff is wanted on the arcsin scale too.
EDGE_THRESHOLD_RAW    <- 2
EDGE_THRESHOLD_ARCSIN <- 0

# Fixed edge-width legend breaks for the raw #Synapses plots.
WEIGHT_BREAKS_RAW <- c(0, 3, 5, 10, 15, 20, 25, 30)

# ---- One plot per resolution/weight variant: circuit minus A19f/A19l,
# Fruchterman-Reingold layout. ----

build_all_combos <- function(mat, row_meta, row_name_col, col_meta, col_name_col,
                              weight_kind, edge_threshold, base_prefix,
                              fig_width, fig_height) {
  build_gephi_network(mat, row_meta, row_name_col, col_meta, col_name_col,
                       weight_kind = weight_kind, edge_threshold = edge_threshold,
                       out_prefix = paste0(base_prefix, "_minusA19_fr"),
                       layout_algo = "fr", exclude_groups = EXCLUDE_A19,
                       fig_width = fig_width, fig_height = fig_height)
}

# 1. side x side (fully resolved both axes) -- raw #Synapses
build_all_combos(mat_side_side_raw, row_meta_side, "Row_label_side", col_meta_side, "Col_label_side",
                  weight_kind = "raw", edge_threshold = EDGE_THRESHOLD_RAW,
                  base_prefix = "network_side_side_raw_forcedirected",
                  fig_width = 34, fig_height = 20)

# 2. side x side -- arcsin(sqrt(proportion))
build_all_combos(mat_side_side_arcsin, row_meta_side, "Row_label_side", col_meta_side, "Col_label_side",
                  weight_kind = "arcsin", edge_threshold = EDGE_THRESHOLD_ARCSIN,
                  base_prefix = "network_side_side_arcsin_forcedirected",
                  fig_width = 34, fig_height = 20)

# 3. seg x seg (both sides averaged) -- raw #Synapses
build_all_combos(mat_seg_seg_raw, row_meta_seg, "Row_label_seg", col_meta_seg, "Col_label_seg",
                  weight_kind = "raw", edge_threshold = EDGE_THRESHOLD_RAW,
                  base_prefix = "network_seg_seg_raw_forcedirected",
                  fig_width = 20, fig_height = 14)

# 4. seg x seg -- arcsin(sqrt(proportion))
build_all_combos(mat_seg_seg_arcsin, row_meta_seg, "Row_label_seg", col_meta_seg, "Col_label_seg",
                  weight_kind = "arcsin", edge_threshold = EDGE_THRESHOLD_ARCSIN,
                  base_prefix = "network_seg_seg_arcsin_forcedirected",
                  fig_width = 20, fig_height = 14)

message("\n==================== DONE ====================")
message("4 force-directed plots written (2 resolutions x 2 weight kinds; circuit minus A19f/A19l; Fruchterman-Reingold), ",
        "each as .pdf/.png/.svg.")
