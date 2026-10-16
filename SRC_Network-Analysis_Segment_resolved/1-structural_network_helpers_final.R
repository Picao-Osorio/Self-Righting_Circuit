## ---------------------------------------------------------------------
## Shared infrastructure for the "Structural properties of the network"
## analysis: section A (centrality)and section B
## 
##
## This file does NOT run any analysis itself -- it is sourced by
## run_structural_network_analysis.R and the three section_*.R files, and
## provides:
##   - circuit-component lookup (new, distinct from the existing numeric
##     hierarchy_rank) and a sensory->motor group display order built from it
##   - the two extra "plain proportion" matrices (mat_*_proportion), which
##     didn't exist before -- only raw #Synapses and arcsin(sqrt(proportion))
##     matrices were built by segment_side_data_prep_3.R. These are derived
##     here from edges_side / long_seg_seg, which already carry the
##     untransformed Proportion column, so segment_side_data_prep_3.R itself
##     does not need to change.
##   - build_structural_graph(): the node/edge/attribute builder used by
##     every section below, with an optional exclude_groups argument for the
##     "full circuit" vs. "minus A19f and A19l" comparison
##   - a small ggplot2 theme and the shared color/shape mappings
##   - generic statistics + plotting helpers reused by every per-node metric
##     in section A: group boxplot+jitter+CLD-letters, paired L/R plot +
##     Wilcoxon, and A-P axis regression plot
##   - an Excel-writing helper
##
## Build_structural_graph() now keeps ONLY SOURCE NEURONS
## (nodes with a row in the connectivity matrix, i.e. sources in the workbook's
## Filtered_downstream sheet) PLUS the LT1/LT2 motoneurons, which are targets
## only but always kept (SOURCE_NODES_ONLY / SOURCE_EXCEPTION_GROUPS below).
## Instances that are purely downstream targets (e.g. A09l (DnB) at T3/A8) are
## dropped together with the edges ending on them. This applies to every
## section that builds graphs here (A, B). Set SOURCE_NODES_ONLY <- FALSE
## to restore the old behaviour.
##
## Packages beyond those 0-segment_side_data_prep_final.R already requires:
## igraph, tidygraph, ggplot2, openxlsx, rstatix, multcompView, broom,
## scales. mclust/aricode/sna/vegan are required later (section B/C only).
## ---------------------------------------------------------------------

required_pkgs_structural <- c(
  "igraph", "tidygraph", "ggplot2", "openxlsx", "rstatix", "multcompView",
  "broom", "scales", "purrr"
)
missing_pkgs_structural <- required_pkgs_structural[
  !vapply(required_pkgs_structural, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs_structural) > 0) install.packages(missing_pkgs_structural)
invisible(lapply(required_pkgs_structural, library, character.only = TRUE))

source("0-segment_side_data_prep_final.R")

## =========================================================================
## 1. Circuit-component lookup (NEW -- distinct from the existing numeric
## hierarchy_rank used for dendrogram/node-position ordering elsewhere).
## =========================================================================


circuit_component_map <- c(
  "mdIV" = "Sensory", "chordotonals" = "Sensory",
  "A09a (Basin-2)" = "1st order IN", "A09g (Basin-3)" = "1st order IN", "A09c (Basin-4)" = "1st order IN",
  "A05q" = "2nd order IN", "A23g" = "2nd order IN",
  "Goro" = "Hub IN", "TePn05" = "Hub IN", "A02o (Wave-1)" = "Hub IN", "A09l (DnB)" = "Hub IN",
  "A27k" = "Prepremotor IN", 
  "A14a" = "Premotor IN", "A01c" = "Premotor IN", "A18j" = "Premotor IN", "A26f" = "Premotor IN",
  "LT1" = "Motoneuron", "LT2" = "Motoneuron"
)

circuit_component_order <- c("Sensory", "1st order IN", "2nd order IN", "Hub IN",
                              "Prepremotor IN", "Premotor IN", "Motoneuron")

get_circuit_component <- function(groups) {
  comp <- circuit_component_map[groups]
  if (any(is.na(comp))) {
    message("No circuit-component assignment for: ", paste(unique(groups[is.na(comp)]), collapse = ", "),
            " -- defaulting to NA (excluded from component-based analyses).")
  }
  unname(comp)
}

# Sensory -> motor display order for every per-group plot in section A,
# built directly from circuit_component_order so the grouping and the
# figure's left-to-right order can never drift apart.
neuron_group_order <- unlist(
  lapply(circuit_component_order, function(comp) names(circuit_component_map)[circuit_component_map == comp]),
  use.names = FALSE
)

## =========================================================================
## 2. Plain-proportion matrices (raw Synapses/Total_Synapses, NOT
## arcsin(sqrt(.))) -- "proportion of synapses" for this
## analysis, distinct from the arcsin-transformed matrices used for the
## heatmap/clustering work. Built from edges_side / long_seg_seg, both of
## which already carry the untransformed Proportion column.
## =========================================================================

mat_side_side_proportion <- edges_side %>%
  select(Row_label_side, Col_label_side, Proportion) %>%
  pivot_wider(names_from = Col_label_side, values_from = Proportion, values_fill = 0) %>%
  column_to_rownames("Row_label_side") %>%
  as.matrix()

mat_seg_seg_proportion <- long_seg_seg %>%
  select(Row_label_seg, Col_label_seg, Proportion_mean) %>%
  pivot_wider(names_from = Col_label_seg, values_from = Proportion_mean, values_fill = 0) %>%
  column_to_rownames("Row_label_seg") %>%
  as.matrix()

## =========================================================================
## 3. Graph builder, shared by every section below.
## =========================================================================
## Mirrors build_topology_graph() in network_topology_analysis.R (including
## its check.names = FALSE fix and explicit NA-metadata guard :
##   - circuit_component / hierarchy_rank_val / segi vertex attributes
##   - an optional exclude_groups argument, for the "full circuit" vs.
##     "minus A19f and A19l" comparison: nodes whose Group is in
##     exclude_groups are dropped BEFORE the graph is built, which also
##     drops every edge touching them (rather than zeroing them out and
##     leaving disconnected/zero-degree nodes behind).

SOURCE_NODES_ONLY       <- TRUE            # TRUE = keep source neurons (+ exception groups) only
SOURCE_EXCEPTION_GROUPS <- c("LT1", "LT2") # motoneurons: targets only in the data, but always kept

build_structural_graph <- function(mat, row_meta, row_name_col, col_meta, col_name_col,
                                    exclude_groups = NULL) {
  meta_to_nodes <- function(meta, name_col) {
    meta <- meta %>% rename(name = !!name_col)
    if (!"Side" %in% names(meta)) meta$Side <- "none"
    meta %>% select(name, Group, Seg, Side)
  }

  # check.names = FALSE: 
  # several canonical Group names contain parentheses/spaces (e.g. "A09a
  # (Basin-2)"), which as.data.frame()'s default check.names = TRUE would
  # silently mangle via make.names(), breaking the left_join() below for
  # every such TARGET node.
  edges <- mat %>%
    as.data.frame(check.names = FALSE) %>%
    rownames_to_column("source") %>%
    pivot_longer(-source, names_to = "target", values_to = "weight") %>%
    filter(!is.na(weight), weight > 0)

  # Source neurons only (+ LT1/LT2): drop target-only nodes and every edge
  # ending on them BEFORE the node table is built (see header note).
  if (SOURCE_NODES_ONLY) {
    exception_names <- col_meta[[col_name_col]][col_meta$Group %in% SOURCE_EXCEPTION_GROUPS]
    source_names <- union(unique(row_meta[[row_name_col]]), exception_names)
    edges <- edges %>% filter(source %in% source_names, target %in% source_names)
  }

  all_nodes_meta <- bind_rows(
    meta_to_nodes(row_meta, row_name_col),
    meta_to_nodes(col_meta, col_name_col)
  ) %>% distinct(name, .keep_all = TRUE)

  node_names <- union(unique(edges$source), unique(edges$target))
  nodes <- tibble(name = node_names) %>%
    left_join(all_nodes_meta, by = "name")

  if (any(is.na(nodes$Group))) {
    stop("Internal error: node(s) with no metadata match: ",
         paste(nodes$name[is.na(nodes$Group)], collapse = ", "))
  }

  if (!is.null(exclude_groups) && length(exclude_groups) > 0) {
    keep_names <- nodes$name[!(nodes$Group %in% exclude_groups)]
    nodes <- nodes %>% filter(name %in% keep_names)
    edges <- edges %>% filter(source %in% keep_names, target %in% keep_names)
  }

  nodes <- nodes %>%
    mutate(
      type               = get_type(Group),
      circuit_component  = get_circuit_component(Group),
      segi               = seg_index(Seg),
      hierarchy_rank_val = get_circuit_rank(Group),
      color              = get_colors_by_group(Group)
    )

  tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
}

## Builds the full set of graph variants this analysis compares: 2 matrix
## resolutions (seg_seg / side_side) x 2 weight kinds (raw / proportion) x
## 2 circuit versions (full / minus_A19) = 8 graphs, returned as a flat
## named list keyed "<resolution>_<weight>_<circuit>".

EXCLUDE_A19 <- c("A19f", "A19l") # neural inhibition of these neurons does not impact SR timing

build_all_structural_graphs <- function() {
  specs <- list(
    seg_seg_raw          = list(mat_seg_seg_raw,          row_meta_seg,  "Row_label_seg",  col_meta_seg,  "Col_label_seg"),
    seg_seg_proportion   = list(mat_seg_seg_proportion,    row_meta_seg,  "Row_label_seg",  col_meta_seg,  "Col_label_seg"),
    side_side_raw        = list(mat_side_side_raw,         row_meta_side, "Row_label_side", col_meta_side, "Col_label_side"),
    side_side_proportion = list(mat_side_side_proportion,  row_meta_side, "Row_label_side", col_meta_side, "Col_label_side")
  )
  out <- list()
  for (nm in names(specs)) {
    s <- specs[[nm]]
    out[[paste0(nm, "_full")]] <- build_structural_graph(s[[1]], s[[2]], s[[3]], s[[4]], s[[5]], exclude_groups = NULL)
    out[[paste0(nm, "_minusA19")]] <- build_structural_graph(s[[1]], s[[2]], s[[3]], s[[4]], s[[5]], exclude_groups = EXCLUDE_A19)
  }
  out
}

## =========================================================================
## 4. Shared aesthetics: 
## =========================================================================
## FLAG (shapes): base ggplot/R pch has no true hexagon marker. Sensory
## stays a triangle (pch 24/ggplot shape 24) and motoneuron a circle (shape
## 21), matching the existing network plots; interneuron uses a diamond
## (shape 23) as the closest built-in stand-in for "rotated hexagon" -- if
## an exact hexagon is wanted, install the 'ggstar' package and swap
## geom_point(shape=...) for ggstar::geom_star(starshape = "hexagon2") in
## the plotting functions below.

shape_for_type <- c(sensory = 24, interneuron = 23, motoneuron = 21)

theme_prism_like <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line         = element_line(colour = "black", linewidth = 0.6),
      axis.ticks        = element_line(colour = "black", linewidth = 0.6),
      axis.text         = element_text(colour = "black", size = 15),
      axis.text.x       = element_text(angle = 45, hjust = 1),
      panel.grid        = element_blank(),
      legend.key        = element_blank(),
      plot.title        = element_text(face = "bold", size = base_size),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold")
    )
}

## =========================================================================
## 5. Group-level statistics: normality -> ANOVA+Bonferroni OR
## Kruskal-Wallis+Dunn (Bonferroni), then compact-letter-display (CLD).
## =========================================================================
## Returns list(test = "anova"/"kruskal", overall_p, pairwise = tibble,
## letters = named character vector keyed by group).

run_group_test_with_letters <- function(df, value_col, group_col) {
  df <- df %>% filter(!is.na(.data[[value_col]]))
  groups <- unique(df[[group_col]])

  # Per-group Shapiro-Wilk (needs >= 3 and <= 5000 obs; groups outside that
  # range are skipped and treated as "not confirmed normal", i.e. push
  # toward the non-parametric branch -- the conservative choice).
  shapiro_ps <- vapply(groups, function(gr) {
    x <- df[[value_col]][df[[group_col]] == gr]
    if (length(x) >= 3 && length(x) <= 5000 && length(unique(x)) > 1) {
      suppressWarnings(stats::shapiro.test(x)$p.value)
    } else {
      NA_real_
    }
  }, numeric(1))

  all_normal <- length(shapiro_ps) > 0 && all(!is.na(shapiro_ps)) && all(shapiro_ps > 0.05)

  if (all_normal) {
    fit <- stats::aov(stats::as.formula(paste0(value_col, " ~ ", group_col)), data = df)
    overall_p <- summary(fit)[[1]][["Pr(>F)"]][1]
    pw <- stats::pairwise.t.test(df[[value_col]], df[[group_col]], p.adjust.method = "bonferroni")
    pmat <- pw$p.value
    pairwise_tbl <- tibble::as_tibble(as.data.frame(as.table(pmat)), .name_repair = "minimal")
    names(pairwise_tbl) <- c("group1", "group2", "p_adj")
    pairwise_tbl <- pairwise_tbl %>% filter(!is.na(p_adj)) %>% mutate(method = "pairwise t-test (Bonferroni)")
    test_used <- "one-way ANOVA + pairwise t-tests (Bonferroni)"
  } else {
    overall_p <- stats::kruskal.test(stats::as.formula(paste0(value_col, " ~ ", group_col)), data = df)$p.value
    dt <- rstatix::dunn_test(df, stats::as.formula(paste0(value_col, " ~ ", group_col)), p.adjust.method = "bonferroni")
    # Reshape rstatix::dunn_test's long output into a square p-value matrix
    # for multcompView::multcompLetters (same shape pairwise.t.test gives).
    all_g <- sort(unique(c(dt$group1, dt$group2)))
    pmat <- matrix(NA_real_, nrow = length(all_g) - 1, ncol = length(all_g) - 1,
                    dimnames = list(all_g[-1], all_g[-length(all_g)]))
    for (i in seq_len(nrow(dt))) {
      g1 <- dt$group1[i]; g2 <- dt$group2[i]
      if (g1 %in% rownames(pmat) && g2 %in% colnames(pmat)) pmat[g1, g2] <- dt$p.adj[i]
      if (g2 %in% rownames(pmat) && g1 %in% colnames(pmat)) pmat[g2, g1] <- dt$p.adj[i]
    }
    # NOTE: rstatix::dunn_test()'s own column is literally named "p.adj"
    # (with a dot, not an underscore) -- selecting "p_adj" doesn't exist
    # and threw "Can't select columns that don't exist" on every single
    # non-normal metric (which, in practice, is most of them), which is
    # why no letters were showing up again. Renamed on the way in instead.
    pairwise_tbl <- dt %>% select(group1, group2, p_adj = p.adj) %>% mutate(method = "Dunn's test (Bonferroni)")
    test_used <- "Kruskal-Wallis + Dunn's test (Bonferroni)"
  }

  # multcompView wants a flat named vector "g1-g2" = p (the SEPARATOR IS
  # HARDCODED to "-" inside multcompLetters -- there is no argument to
  # change it; `compare` is a different parameter entirely: it's the
  # comparison OPERATOR run against `threshold` via get(compare), e.g.
  # "<" or "<=" -- NOT a name separator. Passing "::" there last time made
  # multcompLetters call get("::")(x, threshold), i.e. try to use the
  # namespace-access operator as a function, which is what actually threw
  # "bad namespace name".
  #
  # The real problem is that several of this circuit's canonical group
  # names contain a literal "-" themselves (e.g. "A02o (Wave-1)"), so a
  # pair name like "A02o (Wave-1)-A05q" gets mis-split on the wrong hyphen
  # by the hardcoded "-" separator. Since the separator can't be changed,
  # the fix instead avoids ever handing multcompLetters a group name with
  # a "-" in it: build the pair names from short hyphen-free placeholder
  # codes ("g1", "g2", ...), run multcompLetters on those, then map the
  # resulting letters back to the real group names afterward.
  group_code <- setNames(paste0("grp", seq_along(groups)), as.character(groups))
  flat_p <- c()
  for (r in rownames(pmat)) for (cl in colnames(pmat)) {
    if (!is.na(pmat[r, cl])) flat_p[paste0(group_code[[r]], "-", group_code[[cl]])] <- pmat[r, cl]
  }
  cld <- if (length(flat_p) > 0) {
    cld_coded <- multcompView::multcompLetters(flat_p, threshold = 0.05)$Letters
    code_to_group <- setNames(names(group_code), group_code)
    setNames(unname(cld_coded), code_to_group[names(cld_coded)])
  } else {
    setNames(rep("a", length(groups)), groups)
  }

  list(test = test_used, overall_p = overall_p, letters = cld,
       pairwise = if (exists("pairwise_tbl")) pairwise_tbl else tibble())
}

## =========================================================================
## 6. Plotting helpers (reused identically for every section-A metric --
## degree, strength, betweenness, eigenvector/PageRank, hub/authority,
## closeness, path length, eccentricity).
## =========================================================================

## 6a. Per-group boxplot (transparent, behind) + jitter (in front, colored
## by group / shaped by type) + CLD significance letters above each box.


plot_group_box_jitter <- function(df, value_col, group_col = "Group", type_col = "type",
                                   letters = NULL, y_lab, title = NULL,
                                   group_order = neuron_group_order) {
  df[[group_col]] <- factor(df[[group_col]], levels = intersect(group_order, unique(df[[group_col]])))
  group_colors <- get_colors_by_group(levels(df[[group_col]]))
  names(group_colors) <- levels(df[[group_col]])

  p <- ggplot(df, aes(x = .data[[group_col]], y = .data[[value_col]])) +
    geom_boxplot(aes(fill = .data[[group_col]]), alpha = 0.25, outlier.shape = NA, width = 0.6) +
    geom_jitter(aes(color = .data[[group_col]], shape = .data[[type_col]]),
                width = 0.15, size = 1.8, alpha = 0.85) +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_color_manual(values = group_colors, guide = "none") +
    scale_shape_manual(values = shape_for_type, guide = "none") +
    labs(x = NULL, y = y_lab, title = title) +
    theme_prism_like()

  if (!is.null(letters) && length(letters) > 0) {
    y_max <- df %>% group_by(.data[[group_col]]) %>% summarise(ymax = max(.data[[value_col]], na.rm = TRUE), .groups = "drop")
    y_rng <- diff(range(df[[value_col]], na.rm = TRUE))
    y_max$label <- unname(letters[as.character(y_max[[group_col]])])
    y_max$label[is.na(y_max$label)] <- ""
    p <- p + geom_text(data = y_max, aes(x = .data[[group_col]], y = ymax + 0.06 * y_rng, label = label),
                        inherit.aes = FALSE, fontface = "bold", size = 3.2)
  }
  p
}

## 6b. Side-side only: paired Left/Right plot per neuron group (points
## connected by a grey line per matched segment) + transparent boxplot
## behind + Wilcoxon matched-pairs p-value annotated per group facet.
## Matches Fig.1E/H/K and Fig.2C/F.
##
## `df` must have one row per (Group, Seg, Side) node with columns Group,
## Seg, Side, and value_col. Pairs are matched on (Group, Seg); a Group/Seg
## combination missing either L or R is dropped from the paired test (and
## flagged) but both points are still shown if present.

plot_paired_LR <- function(df, value_col, group_order = neuron_group_order) {
  df <- df %>% filter(Side %in% c("L", "R"))
  wide <- df %>%
    select(Group, Seg, Side, value = all_of(value_col)) %>%
    pivot_wider(names_from = Side, values_from = value)

  incomplete <- wide %>% filter(is.na(L) | is.na(R))
  if (nrow(incomplete) > 0) {
    message(nrow(incomplete), " (Group, Seg) combination(s) missing L or R are excluded from the ",
            "paired Wilcoxon test (point still plotted if present): ",
            paste(paste0(incomplete$Group, "_", incomplete$Seg), collapse = ", "))
  }
  paired <- wide %>% filter(!is.na(L), !is.na(R))

  wilcox_results <- paired %>%
    group_by(Group) %>%
    summarise(
      n = n(),
      p_value = if (n() >= 3 && !all(L == R)) {
        suppressWarnings(stats::wilcox.test(L, R, paired = TRUE)$p.value)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(label = case_when(
      is.na(p_value)  ~ "n/a (n<3)",
      p_value < 0.001 ~ "****",
      p_value < 0.01  ~ "***",
      p_value < 0.05  ~ "**",
      p_value < 0.1   ~ "*",
      TRUE            ~ "n.s."
    ))

  plot_df <- df %>% filter(Group %in% unique(wide$Group))
  plot_df$Group <- factor(plot_df$Group, levels = intersect(group_order, unique(plot_df$Group)))
  group_colors <- get_colors_by_group(levels(plot_df$Group))
  names(group_colors) <- levels(plot_df$Group)

  pair_lines <- paired %>% pivot_longer(c(L, R), names_to = "Side", values_to = "value")
  pair_lines$Group <- factor(pair_lines$Group, levels = levels(plot_df$Group))

  p <- ggplot(plot_df, aes(x = Side, y = .data[[value_col]])) +
    geom_boxplot(aes(fill = Group), alpha = 0.25, outlier.shape = NA, width = 0.5) +
    geom_line(data = pair_lines, aes(x = Side, y = value, group = Seg), color = "grey60", linewidth = 0.3) +
    geom_point(data = pair_lines, aes(x = Side, y = value, color = Group), size = 1.6) +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_color_manual(values = group_colors, guide = "none") +
    facet_wrap(~Group, nrow = 1) +
    labs(x = NULL, y = value_col) +
    theme_prism_like() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  # Annotate each facet with its Wilcoxon result.
  ann <- wilcox_results
  ann$Group <- factor(ann$Group, levels = levels(plot_df$Group))
  y_top <- max(df[[value_col]], na.rm = TRUE)
  p <- p + geom_text(data = ann, aes(x = 1.5, y = y_top * 1.05, label = label),
                      inherit.aes = FALSE, fontface = "bold", size = 3)

  list(plot = p, stats = wilcox_results)
}

## 6c. Along the A-P axis: scatter (per-group color/shape) + per-group lm
## fit line + 95% CI ribbon, plus an overall (all-groups-pooled) fit.
## Slope p-value and R^2 annotated per group, bold when p < 0.05. Matches
## Fig.2D/G and Fig.4D.

plot_AP_regression <- function(df, value_col, group_col = "Group", seg_col = "segi",
                                group_order = neuron_group_order) {
  df <- df %>% filter(!is.na(.data[[value_col]]), !is.na(.data[[seg_col]]))
  df[[group_col]] <- factor(df[[group_col]], levels = intersect(group_order, unique(df[[group_col]])))
  group_colors <- get_colors_by_group(levels(df[[group_col]]))
  names(group_colors) <- levels(df[[group_col]])

  fit_one <- function(sub) {
    if (nrow(sub) < 3 || length(unique(sub[[seg_col]])) < 2) {
      return(tibble(p_value = NA_real_, r_squared = NA_real_, slope = NA_real_))
    }
    m <- stats::lm(stats::as.formula(paste0(value_col, " ~ ", seg_col)), data = sub)
    s <- summary(m)
    tibble(
      p_value   = s$coefficients[2, 4],
      r_squared = s$r.squared,
      slope     = s$coefficients[2, 1]
    )
  }

  per_group_fit <- df %>% group_by(.data[[group_col]]) %>% group_modify(~fit_one(.x)) %>% ungroup()
  overall_fit <- fit_one(df)

  per_group_fit <- per_group_fit %>%
    mutate(label = sprintf("%s: p = %s, R² = %.2f", .data[[group_col]],
                            ifelse(is.na(p_value), "n/a", formatC(p_value, digits = 4, format = "g")),
                            r_squared),
           bold = !is.na(p_value) & p_value < 0.05,
           # Legend-entry version of the same text (used below so the
           # per-group stats sit on the right, inside the legend itself,
           # rather than as a separate text block over the data) --
           # on the right side together with the legend. " *" flags
           # p < 0.05 (legend text can't selectively bold one entry).
           legend_label = sprintf("%s: p = %s, R² = %.2f%s", .data[[group_col]],
                                    ifelse(is.na(p_value), "n/a", formatC(p_value, digits = 4, format = "g")),
                                    r_squared, ifelse(bold, " *", "")))

  # geom_smooth(..., se = TRUE) is already fit PER GROUP here (color/fill are
  # mapped to group_col, and ggplot2 fits one smooth per discrete
  # color/fill group automatically) -- but letting it try on every group
  # regardless of group size is what produced the "qt(): NaNs produced" /
  # "max(ids, na.rm=TRUE): ...-Inf" warning storm: a group with only 2
  # points has residual df = 0, so its 95% CI ribbon is mathematically
  # undefined (NaN), and a group with <2 distinct segment positions has no
  # line to fit at all, leaving an empty/degenerate layer that trips
  # ggplot's internal range calculation. Those groups are exactly the ones
  # fit_one() above already excludes from per_group_fit (n<3 or <2 unique
  # segments); applying the SAME threshold to what's handed to geom_smooth
  # keeps every point visible while only fitting/drawing a line+ribbon
  # where it's statistically meaningful, which removes the warnings instead
  # of just suppressing them.
  smoothable_groups <- per_group_fit[[group_col]][!is.na(per_group_fit$p_value)]
  df_smooth <- df %>% filter(.data[[group_col]] %in% smoothable_groups)

  # Put the per-group p/R² INTO the legend itself (right side, next to the
  # color key), rather than as a separate text block over the data.
  # Done by relabeling the color/fill grouping variable's
  # levels to include each group's stats, so ggplot's own legend renders
  # them; the underlying group identity (color, smooth fit) is unchanged.
  group_levels   <- levels(df[[group_col]])
  legend_lookup  <- setNames(per_group_fit$legend_label, as.character(per_group_fit[[group_col]]))
  legend_levels  <- unname(legend_lookup[group_levels])
  df$group_label <- factor(unname(legend_lookup[as.character(df[[group_col]])]), levels = legend_levels)
  if (nrow(df_smooth) > 0) {
    df_smooth$group_label <- factor(unname(legend_lookup[as.character(df_smooth[[group_col]])]), levels = legend_levels)
  }
  legend_colors <- setNames(unname(group_colors[group_levels]), legend_levels)

  p <- ggplot(df, aes(x = .data[[seg_col]], y = .data[[value_col]], color = group_label)) +
    geom_point(aes(shape = type), alpha = 0.5, size = 1.8)
  if (nrow(df_smooth) > 0) {
    p <- p + geom_smooth(data = df_smooth, aes(x = .data[[seg_col]], y = .data[[value_col]], color = group_label, fill = group_label),
                          method = "lm", formula = y ~ x, se = TRUE, alpha = 0.15, linewidth = 0.8)
  }

  p <- p +
    scale_color_manual(values = legend_colors, name = "Group (slope fit; * = p < 0.05)") +
    scale_fill_manual(values = legend_colors, guide = "none") +
    scale_shape_manual(values = shape_for_type, guide = "none") +
    labs(x = "Segment (axial index)", y = value_col,
         subtitle = sprintf("Overall (pooled across groups): p = %s, R² = %.2f",
                             ifelse(is.na(overall_fit$p_value), "n/a", formatC(overall_fit$p_value, digits = 4, format = "g")),
                             overall_fit$r_squared)) +
    theme_prism_like() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "right",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8, face = "bold"),
          legend.key.height = unit(0.9, "lines"))

  list(plot = p, per_group_fit = per_group_fit, overall_fit = overall_fit)
}

## =========================================================================
## 7. Excel-writing helper -- one workbook per metric, one sheet per graph
## variant, each sheet holding both the raw and normalised values side by
## side plus the node's Group/Seg/Side/type/circuit_component metadata.
## =========================================================================

write_metric_excel <- function(results_list, out_file) {
  wb <- openxlsx::createWorkbook()
  for (nm in names(results_list)) {
    sheet_nm <- substr(nm, 1, 31)  # Excel's 31-char sheet-name limit
    openxlsx::addWorksheet(wb, sheet_nm)
    openxlsx::writeData(wb, sheet_nm, results_list[[nm]])
  }
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
  message("Saved ", out_file, " (", length(results_list), " sheet(s)).")
}

message("structural_network_helpers.R loaded: ", length(circuit_component_map), " groups mapped to ",
        length(circuit_component_order), " circuit components; built ", 2, " proportion matrices; ",
        "build_all_structural_graphs() will build 8 graph variants (2 resolutions x 2 weight kinds x ",
        "2 circuit versions).")
