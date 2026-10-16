## ---------------------------------------------------------------------
## Section A: Centrality properties.
## Sources 1-structural_network_helpers_final.R, builds the graph variants and
## keeps only the "minus A19f/A19l" ones, which are also
## restricted to SOURCE neurons + the LT1/LT2 motoneurons (done in
## build_structural_graph(); target-only nodes such as DnB at T3/A8 are gone).
## LT1/LT2 have no outgoing connections, so they are EXCLUDED from the
## out_degree and out_strength tables and plots (see
## PER_NODE_METRIC_EXCLUDE_GROUPS below). Computes every
## per-node centrality metric, and for each one:
##   - writes raw + normalised values to one Excel workbook (one sheet per
##     graph variant)
##   - draws the three plot types from the brief, generically, for every
##     metric: (i) per-group boxplot+jitter+significance letters,
##     (ii) side_side-only paired L/R plot + Wilcoxon, (iii) A-P axis
##     regression plot
## Graph-level scalars (centralization, reciprocity) and the two
## restricted-vertex-set metrics (path length, eccentricity) are handled
## separately at the end, since they don't fit the generic per-node engine.
##
## Output layout: ./structural_output/A_centrality/<metric>/
## ---------------------------------------------------------------------

source("1-structural_network_helpers_final.R")

OUT_DIR_A <- "structural_output/A_centrality"
dir.create(OUT_DIR_A, recursive = TRUE, showWarnings = FALSE)

message("Building the graph variants...")
graphs <- build_all_structural_graphs()
# Keep only the "minus A19f/A19l" circuit version.
full_variant_names <- names(graphs)
graphs <- graphs[grepl("minusA19", names(graphs))]
if (length(graphs) == 0) {
  stop("Section A: no 'minusA19' variant found among ", paste(full_variant_names, collapse = ", "))
}
message("Section A: keeping only minusA19 variants (", length(graphs), " of ", length(full_variant_names), ").")
for (nm in names(graphs)) {
  message("  ", nm, ": ", igraph::gorder(graphs[[nm]]), " nodes, ", igraph::gsize(graphs[[nm]]), " edges.")
}

## =========================================================================
## 1A/2A/3A/4A/5A: generic per-node metric engine.
## =========================================================================
## Each entry: function(g) -> list(raw = named numeric vector, norm = named
## numeric vector), both named by node (V(g)$name).

PER_NODE_METRICS <- list(
  ## ---- 1A: degree and strength --------------------------------------
  in_degree = function(g) {
    r <- igraph::degree(g, mode = "in")
    list(raw = r, norm = r / (igraph::gorder(g) - 1))
  },
  out_degree = function(g) {
    r <- igraph::degree(g, mode = "out")
    list(raw = r, norm = r / (igraph::gorder(g) - 1))
  },
  in_strength = function(g) {
    r <- igraph::strength(g, mode = "in", weights = igraph::E(g)$weight)
    list(raw = r, norm = if (max(r) > 0) r / max(r) else r)
  },
  out_strength = function(g) {
    r <- igraph::strength(g, mode = "out", weights = igraph::E(g)$weight)
    list(raw = r, norm = if (max(r) > 0) r / max(r) else r)
  },
  ## ---- 2A: betweenness ------------------------------------------------
  ## Path "distance" = 1/weight (a strong connection = a short/cheap path);
  ## using raw weight directly would backwards-treat well-connected pairs
  ## as far apart. The unweighted (topological) version is included too,
  ## as explicitly requested.
  betweenness_weighted = function(g) {
    list(
      raw  = igraph::betweenness(g, directed = TRUE, weights = 1 / igraph::E(g)$weight, normalized = FALSE),
      norm = igraph::betweenness(g, directed = TRUE, weights = 1 / igraph::E(g)$weight, normalized = TRUE)
    )
  },
  betweenness_unweighted = function(g) {
    list(
      raw  = igraph::betweenness(g, directed = TRUE, weights = NA, normalized = FALSE),
      norm = igraph::betweenness(g, directed = TRUE, weights = NA, normalized = TRUE)
    )
  },
  ## ---- 3A: eigenvector centrality / PageRank --------------------------
  eigenvector = function(g) {
    list(
      raw  = igraph::eigen_centrality(g, directed = TRUE, weights = igraph::E(g)$weight, scale = FALSE)$vector,
      norm = igraph::eigen_centrality(g, directed = TRUE, weights = igraph::E(g)$weight, scale = TRUE)$vector
    )
  },
  pagerank = function(g) {
    r <- igraph::page_rank(g, directed = TRUE, weights = igraph::E(g)$weight, damping = 0.85)$vector
    list(raw = r, norm = if (max(r) > 0) r / max(r) else r)  # raw sums to 1 by construction; norm rescales to max=1 for comparability with the other metrics here
  },
  ## ---- 4A: hub / authority (Kleinberg HITS) ---------------------------
  hub_score = function(g) {
    list(
      raw  = igraph::hub_score(g, weights = igraph::E(g)$weight, scale = FALSE)$vector,
      norm = igraph::hub_score(g, weights = igraph::E(g)$weight, scale = TRUE)$vector
    )
  },
  authority_score = function(g) {
    list(
      raw  = igraph::authority_score(g, weights = igraph::E(g)$weight, scale = FALSE)$vector,
      norm = igraph::authority_score(g, weights = igraph::E(g)$weight, scale = TRUE)$vector
    )
  },
  ## ---- 5A: closeness ---------------------------------------------------
  ## FLAG: this is a near-feedforward, largely acyclic directed circuit
  ## (sensory -> ... -> motor), so most node pairs are NOT mutually
  ## reachable in both directions. mode = "out" (how quickly a node can
  ## REACH the rest of the network) is used here rather than igraph's
  ## default "all", because "all" silently treats the directed graph as
  ## undirected for reachability, which defeats the point of a directed
  ## circuit. Expect many small/near-equal values for downstream
  ## (motoneuron-adjacent) nodes, which can reach almost nothing -- this
  ## metric is intrinsically a poor fit for this network's topology (see
  ## the written summary for a fuller discussion).
  closeness_out = function(g) {
    list(
      raw  = igraph::closeness(g, mode = "out", weights = 1 / igraph::E(g)$weight, normalized = FALSE),
      norm = igraph::closeness(g, mode = "out", weights = 1 / igraph::E(g)$weight, normalized = TRUE)
    )
  }
)

## Per-metric groups to leave out of the per-node tables/plots: the LT1/LT2
## motoneurons have no outgoing connections, so their out-degree/out-strength
## is structurally 0 and is not an informative measurement.
## Normalisation is still computed on the full graph before the rows are dropped.
PER_NODE_METRIC_EXCLUDE_GROUPS <- list(
  out_degree   = c("LT1", "LT2"),
  out_strength = c("LT1", "LT2")
)

## Neuron categories for the A-P regression plots. A19f and
## A19l are not listed because they are excluded from this analysis.
AP_CATEGORIES <- list(
  `1_sensory_basins` = c("mdIV", "chordotonals", "A09a (Basin-2)", "A09g (Basin-3)", "A09c (Basin-4)"),
  `2_2ndIN_hub`      = c("A05q", "A23g", "Goro", "TePn05"),
  `3_hub2`           = c("A02o (Wave-1)", "A09l (DnB)", "A27k"),
  `4_premotor`       = c("A14a", "A01c", "A18j", "A26f"),
  `5_motoneurons`    = c("LT1", "LT2")
)
## Neuron colours used ONLY in the A-P regression plots;
## all other plots keep the project-wide neuron_color_map. A19f/A19l are not
## listed because they are excluded from this analysis.
neuron_color_map_APregression <- c(
  "LT1"             = "#FF9300",
  "LT2"             = "#0000FF",
  "A26f"            = "#996733",
  "A18j"            = "#03ff00",
  "A01c"            = "#800080",
  "A14a"            = "#ff6666",
  "A09l (DnB)"      = "#ff2501",
  "A27k"            = "#942093",
  "A02o (Wave-1)"   = "#800000",
  "TePn05"          = "#400040",
  "Goro"            = "#2c3f50",
  "A23g"            = "#ffe97f",
  "A05q"            = "#ffd479",
  "A09c (Basin-4)"  = "#008080",
  "A09g (Basin-3)"  = "#008f01",
  "A09a (Basin-2)"  = "#7e8024",
  "chordotonals"    = "#02ffff",
  "mdIV"            = "#013f80"
)

AP_CATEGORY_LABELS <- list(
  `1_sensory_basins` = "sensory + basins", `2_2ndIN_hub` = "2nd-order IN + hub",
  `3_hub2` = "hub 2 (Wave-1, DnB, A27k)", `4_premotor` = "premotor", `5_motoneurons` = "motoneurons"
)

nodes_df <- function(g) {
  g %>% tidygraph::activate(nodes) %>% as_tibble() %>%
    select(name, Group, Seg, Side, type, circuit_component, segi, hierarchy_rank_val)
}

run_per_node_metric <- function(metric_name, compute_fn) {
  message("Computing ", metric_name, " across all ", length(graphs), " variants...")
  excl_groups <- PER_NODE_METRIC_EXCLUDE_GROUPS[[metric_name]]
  out_dir <- file.path(OUT_DIR_A, metric_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  excel_sheets <- list()
  combined <- list()

  for (variant in names(graphs)) {
    g <- graphs[[variant]]
    res <- compute_fn(g)
    nd <- nodes_df(g)
    nd[[paste0(metric_name, "_raw")]]  <- unname(res$raw[nd$name])
    nd[[paste0(metric_name, "_norm")]] <- unname(res$norm[nd$name])
    if (!is.null(excl_groups)) nd <- nd %>% filter(!(Group %in% excl_groups))
    excel_sheets[[variant]] <- nd
    nd$variant <- variant
    combined[[variant]] <- nd
  }
  write_metric_excel(excel_sheets, file.path(out_dir, paste0(metric_name, ".xlsx")))

  all_df <- bind_rows(combined)
  raw_col  <- paste0(metric_name, "_raw")
  norm_col <- paste0(metric_name, "_norm")

  # (i) per-group boxplot+jitter+letters -- one per variant, for both raw
  # and normalised values.
  for (variant in names(graphs)) {
    df_v <- all_df %>% filter(variant == !!variant)
    for (value_col in c(raw_col, norm_col)) {
      if (length(unique(df_v$Group)) < 2 || sum(!is.na(df_v[[value_col]])) < 2) next
      gt <- tryCatch(run_group_test_with_letters(df_v, value_col, "Group"),
                      error = function(e) {
                        message("Group significance test FAILED for ", metric_name, "/", variant, "/", value_col,
                                " -- no letters will be shown on this boxplot. Reason: ", conditionMessage(e))
                        NULL
                      })
      p <- plot_group_box_jitter(df_v, value_col, letters = if (!is.null(gt)) gt$letters else NULL,
                                  y_lab = value_col,
                                  title = paste0(metric_name, " -- ", variant,
                                                 if (!is.null(gt)) paste0(" (", gt$test, ")") else ""))
      ggsave(file.path(out_dir, paste0(metric_name, "_", variant, "_", value_col, "_bygroup.pdf")),
             p, width = 10, height = 5)

      # Collect the test/pairwise/letters results so they end up in an
      # Excel workbook, not just baked into the plot title/annotations --
      # these three lists are declared at the top level (just before the
      # "for (metric_name in names(PER_NODE_METRICS))" loop below) and
      # appended to here via <<- since run_per_node_metric() is defined in
      # that same top-level environment.
      if (!is.null(gt)) {
        group_stats_overall[[length(group_stats_overall) + 1]] <<- tibble(
          metric = metric_name, variant = variant, value_col = value_col,
          test = gt$test, overall_p = gt$overall_p
        )
        if (!is.null(gt$pairwise) && nrow(gt$pairwise) > 0) {
          group_stats_pairwise[[length(group_stats_pairwise) + 1]] <<- gt$pairwise %>%
            mutate(metric = metric_name, variant = variant, value_col = value_col, .before = 1)
        }
        if (length(gt$letters) > 0) {
          group_stats_letters[[length(group_stats_letters) + 1]] <<- tibble(
            metric = metric_name, variant = variant, value_col = value_col,
            Group = names(gt$letters), letter = unname(gt$letters)
          )
        }
      }
    }
  }

  # (ii) side_side-only paired L/R plot -- raw and normalised.
  for (variant in grep("^side_side", names(graphs), value = TRUE)) {
    df_v <- all_df %>% filter(variant == !!variant)
    for (value_col in c(raw_col, norm_col)) {
      if (sum(!is.na(df_v[[value_col]])) < 2) next
      res_lr <- tryCatch(plot_paired_LR(df_v, value_col), error = function(e) NULL)
      if (is.null(res_lr)) next
      ggsave(file.path(out_dir, paste0(metric_name, "_", variant, "_", value_col, "_LRpaired.pdf")),
             res_lr$plot, width = 14, height = 4)
      write.csv(res_lr$stats, file.path(out_dir, paste0(metric_name, "_", variant, "_", value_col, "_LRpaired_wilcoxon.csv")),
                row.names = FALSE)
    }
  }

  # (iii) A-P axis regression -- one plot per variant x value column x
  # neuron CATEGORY (split into 5 categories instead of one
  # plot with every group). Categories are defined by AP_CATEGORIES below; the
  # slope fits and the "overall (pooled)" fit in each plot are computed within
  # that category only. Categories with fewer than 3 data points (e.g. the
  # motoneurons for out_degree/out_strength, which are excluded there) are skipped.
  for (variant in names(graphs)) {
    df_v <- all_df %>% filter(variant == !!variant)
    for (value_col in c(raw_col, norm_col)) {
      for (cat_name in names(AP_CATEGORIES)) {
        df_c <- df_v %>% filter(Group %in% AP_CATEGORIES[[cat_name]])
        if (sum(!is.na(df_c[[value_col]])) < 3) next
        res_ap <- tryCatch(plot_AP_regression(df_c, value_col), error = function(e) NULL)
        if (is.null(res_ap)) next
        # x axis: real segment names in axial order (T2, T3, A1...A8), one
        # labelled tick per segment -- not just the numeric axial index
        # The full T2-A8 range is always shown, even where a
        # category has no data in some segments. Any segment present in the data outside T2-A8
        # is added so no point is left without a label.
        seg_breaks <- sort(unique(c(2:11, df_c$segi[!is.na(df_c$segi)])))
        seg_labels <- ifelse(seg_breaks <= length(segment_levels), segment_levels[pmin(seg_breaks, length(segment_levels))], "?")
        # Recolour with the A-P regression palette. The helper's legend entries
        # carry each group's stats, so map legend_label -> Group -> colour.
        pg <- res_ap$per_group_fit
        ap_cols <- setNames(unname(neuron_color_map_APregression[as.character(pg$Group)]), pg$legend_label)
        ap_cols[is.na(ap_cols)] <- "grey60"
        res_ap$plot <- suppressMessages(res_ap$plot +
          scale_color_manual(values = ap_cols, name = "Group (slope fit; * = p < 0.05)") +
          scale_fill_manual(values = ap_cols, guide = "none"))
        res_ap$plot <- res_ap$plot +
          scale_x_continuous(breaks = seg_breaks, labels = seg_labels,
                             limits = c(min(seg_breaks) - 0.5, max(seg_breaks) + 0.5)) +
          labs(x = "Segment",
               title = paste0(metric_name, " -- ", variant, " -- ", AP_CATEGORY_LABELS[[cat_name]]))
        base_nm <- paste0(metric_name, "_", variant, "_", value_col, "_APregression_", cat_name)
        ggsave(file.path(out_dir, paste0(base_nm, ".pdf")), res_ap$plot, width = 8, height = 6)
        write.csv(res_ap$per_group_fit, file.path(out_dir, paste0(base_nm, "_fits.csv")), row.names = FALSE)
      }
    }
  }

  message("Done: ", metric_name, " (", length(graphs), " variants, outputs in ", out_dir, ")")
  invisible(all_df)
}

## Accumulators for the 1A-5A group-significance results (Shapiro ->
## ANOVA+Bonferroni or Kruskal-Wallis+Dunn(Bonferroni) -> CLD letters),
## filled in by run_per_node_metric() via <<- and written to a single
## Excel workbook below, so the numbers behind every boxplot's letters are
## available to check/cite, not just baked into the plot image.
group_stats_overall  <- list()
group_stats_pairwise <- list()
group_stats_letters  <- list()

for (metric_name in names(PER_NODE_METRICS)) {
  run_per_node_metric(metric_name, PER_NODE_METRICS[[metric_name]])
}

group_stats_overall_tbl  <- bind_rows(group_stats_overall)
group_stats_pairwise_tbl <- bind_rows(group_stats_pairwise)
group_stats_letters_tbl  <- bind_rows(group_stats_letters)

wb_stats <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_stats, "overall_test"); openxlsx::writeData(wb_stats, "overall_test", group_stats_overall_tbl)
openxlsx::addWorksheet(wb_stats, "pairwise_comparisons"); openxlsx::writeData(wb_stats, "pairwise_comparisons", group_stats_pairwise_tbl)
openxlsx::addWorksheet(wb_stats, "cld_letters"); openxlsx::writeData(wb_stats, "cld_letters", group_stats_letters_tbl)
openxlsx::saveWorkbook(wb_stats, file.path(OUT_DIR_A, "group_significance_tests_1A-5A.xlsx"), overwrite = TRUE)
message("Saved ", file.path(OUT_DIR_A, "group_significance_tests_1A-5A.xlsx"),
        " (", nrow(group_stats_overall_tbl), " metric/variant/value_col test(s), ",
        nrow(group_stats_pairwise_tbl), " pairwise comparison(s), ",
        nrow(group_stats_letters_tbl), " group letter assignment(s)).")



message("Section A (centrality) complete. Outputs in ", OUT_DIR_A)

