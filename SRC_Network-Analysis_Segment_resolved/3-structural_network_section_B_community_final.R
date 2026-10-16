## ---------------------------------------------------------------------
## Section B: Clustering and community properties.

## the graphs analysed here contain only SOURCE neurons
## plus the LT1/LT2 motoneurons (targets only, always kept), and exclude
## A19f/A19l -- the source-only filter is applied in build_structural_graph()
## (1-structural_network_helpers_v6.R), the A19 exclusion by the minusA19
## filter just below. Target-only instances (e.g. DnB at T3/A8) are not nodes.
##
## 1B: three community-detection methods (Louvain, Leiden, walktrap), number of communities, membership,
##     modularity Q.
## 2B: global + local clustering coefficient (transitivity), local values
##     summarised by each method's communities.
##
## Output layout: ./structural_output/B_community/
## ---------------------------------------------------------------------

source("1-structural_network_helpers_final.R")
if (!exists("graphs")) graphs <- build_all_structural_graphs()

# ---- Keep only the "minus A19f/A19l" circuit version:
# drop every variant whose name doesn't mark it as the minusA19 version.
# Done here, once, right after the graphs are built/loaded, so every
# downstream step (1B, 2B, both naturally iterate over names(graphs)/
# unique(community_membership_tbl$variant)) automatically only ever sees
# the minusA19 variants without needing separate filtering logic later. ----
full_variant_names <- names(graphs)
graphs <- graphs[grepl("minusA19", names(graphs))]
if (length(graphs) == 0) {
  stop("Section B: no 'minusA19' variant found among ", paste(full_variant_names, collapse = ", "),
       " -- check build_all_structural_graphs()'s naming convention.")
}
message("Section B: keeping only minusA19 variants (", length(graphs), " of ", length(full_variant_names),
        "): ", paste(names(graphs), collapse = ", "))

OUT_DIR_B <- "structural_output/B_community"
dir.create(OUT_DIR_B, recursive = TRUE, showWarnings = FALSE)

nodes_df <- function(g) {
  g %>% tidygraph::activate(nodes) %>% as_tibble() %>%
    select(name, Group, Seg, Side, type, circuit_component, segi, hierarchy_rank_val)
}

## =========================================================================
## 1B: Community detection -- three methods (Louvain, Leiden, walktrap).
## 
## FLAG: cluster_louvain() and cluster_walktrap() are UNDIRECTED-only in
## igraph, so both are run on an undirected collapse of each graph (parallel
## edges summed, i.e. a mutual A->B and B->A connection becomes one edge with
## their combined weight) -- direction is therefore NOT used by either
## method, and modularity Q is evaluated on that same undirected collapse.
## cluster_walktrap()'s steps: left at igraph's default of 4, which is a
## reasonable, widely-used default for a network this size (~100-300
## nodes); there's no strong a priori reason from this circuit's structure
## to prefer a different value, but it's worth trying steps = 2 and 8 as a
## sensitivity check if the walktrap communities look unstable.
## Reproducibility: cluster_louvain() is stochastic and was never seeded
## (its result depended on how many unseeded random
## draws preceded it in the loop). It now runs from 10 fixed seeds, keeping
## the highest-modularity partition, independent of anything run earlier.
## cluster_walktrap() is deterministic and needs no seed.
## =========================================================================

message("Running community detection (1B) across all minusA19 variants...")

COMMUNITY_SEEDS <- 1:10  # fixed, independent of anything run earlier in the script

run_community_methods <- function(g) {
  g_undir <- igraph::as.undirected(g, mode = "collapse", edge.attr.comb = "sum")

  louvain <- NULL
  for (s in COMMUNITY_SEEDS) {
    set.seed(s)
    cand <- igraph::cluster_louvain(g_undir, weights = igraph::E(g_undir)$weight)
    if (is.null(louvain) || igraph::modularity(cand) > igraph::modularity(louvain)) louvain <- cand
  }

  # Leiden: Louvain's successor (guarantees well-connected
  # communities). Set to objective_function = "modularity" so its Q is
  # comparable with Louvain/walktrap (igraph's default is CPM, which needs a
  # hand-tuned resolution). n_iterations = -1 runs until no further
  # improvement. Same seeds / best-of-N-by-modularity rule as Louvain.
  # The resolution argument is named `resolution` in igraph >= 2.0 and
  # `resolution_parameter` in 1.x; left at its default (1) either way.
  leiden <- NULL
  for (s in COMMUNITY_SEEDS) {
    set.seed(s)
    cand <- igraph::cluster_leiden(g_undir, objective_function = "modularity",
                                   weights = igraph::E(g_undir)$weight,
                                   n_iterations = -1)
    if (is.null(leiden) || igraph::modularity(g_undir, cand$membership, weights = igraph::E(g_undir)$weight) >
        igraph::modularity(g_undir, leiden$membership, weights = igraph::E(g_undir)$weight)) leiden <- cand
  }

  walktrap <- igraph::cluster_walktrap(g_undir, weights = igraph::E(g_undir)$weight, steps = 4)  # deterministic, no seed needed


  list(louvain = louvain, leiden = leiden, walktrap = walktrap)
}

community_summary <- list()
community_membership <- list()

for (variant in names(graphs)) {
  g <- graphs[[variant]]
  methods <- run_community_methods(g)
  nd <- nodes_df(g)

  for (method_name in names(methods)) {
    comm <- methods[[method_name]]
    mem <- igraph::membership(comm)
    # Safe, explicit lookup instead of mem[name] / mem[V(g)$name]: membership()
    # is named by vertex name, but indexing a named vector by a character
    # vector silently returns NA for any name that doesn't match exactly
    # (order mismatch, stray whitespace, or a vertex dropped/renamed
    # upstream), and that NA then reaches aricode::NMI() downstream (2C),
    # which -- unlike mclust::adjustedRandIndex() -- refuses NA input. A
    # join makes any such mismatch explicit here instead of propagating.
    mem_tbl <- tibble(name = names(mem), community = as.integer(mem))

    # modularity(): a single Q value for the WHOLE GRAPH under this one
    # partition (not a per-community score).
    # All methods (louvain/leiden/walktrap) were computed on the undirected
    # collapse, so Q is evaluated on that same collapse, for consistency
    # with what each method optimised.
    g_for_mod <- igraph::as.undirected(g, mode = "collapse", edge.attr.comb = "sum")
    mod_join <- tibble(name = igraph::V(g_for_mod)$name) %>% left_join(mem_tbl, by = "name")
    if (anyNA(mod_join$community)) {
      stop("Section B (1B): ", sum(is.na(mod_join$community)),
           " node(s) in variant '", variant, "', method '", method_name,
           "' got no community assignment when matched by name for modularity(). Missing for: ",
           paste(mod_join$name[is.na(mod_join$community)], collapse = ", "))
    }
    Q <- igraph::modularity(g_for_mod, mod_join$community, weights = igraph::E(g_for_mod)$weight)

    community_summary[[length(community_summary) + 1]] <- tibble(
      variant = variant, method = method_name,
      n_communities = length(unique(mem)), modularity_Q = Q
    )

    nd_join <- nd %>% left_join(mem_tbl, by = "name")
    if (anyNA(nd_join$community)) {
      stop("Section B (1B): ", sum(is.na(nd_join$community)),
           " node(s) in variant '", variant, "', method '", method_name,
           "' got no community assignment when matched by name for membership table. Missing for: ",
           paste(nd_join$name[is.na(nd_join$community)], collapse = ", "))
    }
    community_membership[[length(community_membership) + 1]] <- nd_join %>%
      mutate(variant = variant, method = method_name)
  }
  message("  ", variant, ": done (louvain/leiden/walktrap).")
}

community_summary_tbl <- bind_rows(community_summary)
community_membership_tbl <- bind_rows(community_membership)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "summary"); openxlsx::writeData(wb, "summary", community_summary_tbl)
openxlsx::addWorksheet(wb, "membership"); openxlsx::writeData(wb, "membership", community_membership_tbl)
openxlsx::saveWorkbook(wb, file.path(OUT_DIR_B, "community_detection.xlsx"), overwrite = TRUE)

## Plot 1: modularity Q per method x variant (graph-level comparison bar chart).
p_mod <- ggplot(community_summary_tbl, aes(x = method, y = modularity_Q, fill = method)) +
  geom_col(color = "black") +
  facet_wrap(~variant) +
  labs(x = NULL, y = "Modularity Q", title = "Community-detection modularity by method and variant") +
  theme_prism_like() + theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR_B, "modularity_by_method_variant.pdf"), p_mod, width = 10, height = 7)

## =========================================================================
## Plot 2 (replaces the old nodes x method combined heatmap):
## a "community grid" -- Group on the y-axis (ordered by hierarchy_rank_val,
## i.e. the same sensory->motor ordering used elsewhere in this project),
## segment on the x-axis,
## tiles colored by detected community, ONE PLOT PER DETECTION METHOD (and
## per variant), rather than cramming every method into one plot.
##
## Why this over other options: the earlier heatmap put every node on the
## y-axis (ordered sensory->motor, with segment folded into that same
## axis) and methods side by side as columns -- useful for comparing
## methods against each other, but it mixes segment and group into a
## single dense axis, so it's hard to read "does this community track a
## segmental (rostro-caudal) pattern, a neuron-identity pattern, or both"
## at a glance. A Group x Segment GRID answers exactly that: a community
## that lights up as a clean horizontal band = a neuron-identity-driven
## community (same group, many segments); a clean vertical band = a
## segment-driven community (many groups, one segment); anything blockier
## is a genuine group x segment interaction. This is the standard
## "clustered/annotated heatmap" idiom for exactly this kind of
## bipartite-looking grouping question, just drawn directly with
## ggplot2::geom_tile() rather than a general-purpose clustering-heatmap
## package (pheatmap/ComplexHeatmap), since the row/column ORDER here is
## already fully determined by known biology (hierarchy_rank_val, segment
## index) rather than something to be discovered by the plotting package's
## own dendrogram clustering.
##
## Side handling: for the side x side-resolved variants there are two
## nodes per (Group, Seg) cell (L and R); those are drawn as two side-by-
## side panels (facet_wrap(~Side)) of the same Group x Segment grid rather
## than merged, so the two sides' community assignments stay fully visible
## and comparable. The seg x seg-resolved variants (Side == "none") get a
## single panel, same as before.
## =========================================================================

message("Building per-method community grids (Plot 2)...")

# T1, T3 and A9 are intentionally dropped here 
# ("segments... from T2,A1-A8") -- flagged loudly (not silently) in case
# this was meant to be the full T1-A9 range, since it visibly narrows the
# plotted circuit relative to every other plot in this project.
GRID_SEGS <- c("T2", paste0("A", 1:8))

# A qualitative palette wide enough for a couple dozen communities without
# repeating adjacent hues (RColorBrewer's 12-color "Paired" extended by
# interpolation covers the common case; genuinely huge community counts
# will still cycle, which is an inherent limit of categorical color, not a
# bug here).
community_palette <- function(n) {
  base <- RColorBrewer::brewer.pal(min(max(n, 3), 12), "Paired")
  if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
}

plot_community_grid <- function(df_vm, method_name, variant_name) {
  df_vm <- df_vm %>% filter(Seg %in% GRID_SEGS)
  dropped <- setdiff(unique(df_vm$Seg[df_vm$Seg %in% GRID_SEGS == FALSE]), character(0))  # defensive; always empty after filter
  if (nrow(df_vm) == 0) {
    warning("plot_community_grid: no nodes left for ", method_name, "/", variant_name,
            " after restricting to segments ", paste(GRID_SEGS, collapse = ", "), " -- skipping.")
    return(invisible(NULL))
  }

  # One hierarchy_rank_val per Group (mean, in case of any residual
  # per-node/per-segment variation) to order the y-axis; NA ranks sort
  # last with a loud warning rather than silently dropping those groups.
  group_rank <- df_vm %>%
    group_by(Group) %>%
    summarise(rank_val = mean(hierarchy_rank_val, na.rm = TRUE), .groups = "drop")
  if (anyNA(group_rank$rank_val)) {
    warning("plot_community_grid (", method_name, "/", variant_name, "): ",
            sum(is.na(group_rank$rank_val)), " group(s) have no hierarchy_rank_val -- plotted last on the y-axis: ",
            paste(group_rank$Group[is.na(group_rank$rank_val)], collapse = ", "))
    group_rank$rank_val[is.na(group_rank$rank_val)] <- Inf
  }
  group_levels <- group_rank %>% arrange(rank_val) %>% pull(Group)
  df_vm$Group <- factor(df_vm$Group, levels = rev(group_levels))  # rev(): lowest rank at the TOP of the y-axis
  df_vm$Seg   <- factor(df_vm$Seg, levels = GRID_SEGS)

  n_comm <- dplyr::n_distinct(df_vm$community)
  pal <- setNames(community_palette(n_comm), sort(unique(df_vm$community)))

  p <- ggplot(df_vm, aes(x = Seg, y = Group, fill = factor(community))) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = pal, name = "Community") +
    labs(x = NULL, y = NULL,
         title = paste0("Community grid -- ", method_name, " (", variant_name, ")"),
         subtitle = "Rows ordered by hierarchy_rank_val (sensory -> motor, top to bottom); columns T2, A1-A8") +
    theme_prism_like() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = 7))

  n_sides <- dplyr::n_distinct(df_vm$Side)
  if (n_sides > 1) p <- p + facet_wrap(~Side)

  n_groups <- dplyr::n_distinct(df_vm$Group)
  height <- max(4, 0.28 * n_groups)
  width  <- if (n_sides > 1) 10 else 6

  out_path <- file.path(OUT_DIR_B, paste0("community_grid_", method_name, "_", variant_name, ".pdf"))
  ggsave(out_path, p, width = width, height = height)
  message("  Saved ", out_path)
}

for (variant in unique(community_membership_tbl$variant)) {
  for (method_name in unique(community_membership_tbl$method)) {
    df_vm <- community_membership_tbl %>% filter(variant == !!variant, method == !!method_name)
    plot_community_grid(df_vm, method_name, variant)
  }
}

message("1B (community detection) done.")


####################################
###################################
## Section D: additional structural analyses flagged as worth considering

1D. Null-model significance test for community-detection modularity Q
##       (is Q actually higher than chance, for a degree-and-weight-
##       sequence-preserving random graph?)
## 1D uses a degree-preserving random-rewiring null (igraph::rewire() with
## keeping_degseq()) with the observed edge weights randomly shuffled over the
## rewired edges (rewire() itself drops edge attributes -- see
## rewire_weighted() in 1D). Node degrees are preserved exactly, strengths are
## not. 4D uses the same rewiring but is degree-based, so weights play no role. n_perm/n_null below are kept modest (50-100) to keep runtime
## reasonable; raise them for final numbers (see comments at each call).
##
## Sources 1-structural_network_helpers_final.R (reuses `graphs` if already
## built in this session). 1D recomputes its own community detection with the
## same three methods as section B (louvain, leiden, walktrap).
##
## Output layout: ./structural_output/D_additional/
## ---------------------------------------------------------------------

source("1-structural_network_helpers_final.R")
if (!exists("graphs")) graphs <- build_all_structural_graphs()

# same circuit as sections A/B -- only the "minus A19f/A19l"
# variants, and (via build_structural_graph() in the helpers) only SOURCE
# neurons + the LT1/LT2 motoneurons. This now applies to 2D-4D as well as 1D.
full_variant_names <- names(graphs)
graphs <- graphs[grepl("minusA19", names(graphs))]
if (length(graphs) == 0) stop("Section D: no 'minusA19' variant found among ", paste(full_variant_names, collapse = ", "))
message("Section D: keeping only minusA19 variants (", length(graphs), " of ", length(full_variant_names), ").")

OUT_DIR_D <- "structural_output/D_additional"
dir.create(OUT_DIR_D, recursive = TRUE, showWarnings = FALSE)

nodes_df <- function(g) {
  g %>% tidygraph::activate(nodes) %>% as_tibble() %>%
    select(name, Group, Seg, Side, type, circuit_component, segi, hierarchy_rank_val)
}

## =========================================================================
## 1D: Null-model modularity test.
## Aligned with section B: the three community-detection methods
## are Louvain, Leiden and walktrap (edge_betweenness and directed infomap
## were removed from B -- see the notes in the B script), and the test is run
## on the minusA19(f+l) variants only, as in B.
##
## For each method, is the observed modularity Q higher than for a
## degree-/weight-sequence-preserving random rewiring of the same graph?
## Empirical p-value = (# null Q >= observed Q + 1) / (n_perm + 1).
##
## Like B, all methods run on the UNDIRECTED collapse (parallel edges summed)
## and Q is evaluated on that same collapse. Louvain/Leiden are stochastic:
## each run below uses the SAME best-of-COMMUNITY_SEEDS (highest modularity)
## rule as section B, for the observed graph AND for every null graph, so the
## observed Q is not advantaged over the null by getting extra attempts.
## Every call sets its own seeds, so results are reproducible and do not
## depend on what ran earlier in the session. Leiden uses
## objective_function = "modularity" (igraph's default, CPM, is not
## comparable with Q).
## =========================================================================

message("1D: null-model modularity test (louvain / leiden / walktrap; minusA19 variants)...")

COMMUNITY_SEEDS <- 1:10  # same as section B

best_of_seeds <- function(run_one) {
  best <- NULL; best_q <- -Inf
  for (s in COMMUNITY_SEEDS) {
    set.seed(s)
    cand <- run_one()
    q <- cand$q
    if (q > best_q) { best <- cand$memb; best_q <- q }
  }
  list(membership = best, Q = best_q)
}

# gg: undirected weighted graph. Returns list(membership, Q).
run_community_method <- function(gg, method_name) {
  w <- igraph::E(gg)$weight
  switch(method_name,
    louvain = best_of_seeds(function() {
      cm <- igraph::cluster_louvain(gg, weights = w)
      list(memb = igraph::membership(cm), q = igraph::modularity(gg, igraph::membership(cm), weights = w))
    }),
    leiden = best_of_seeds(function() {
      cm <- igraph::cluster_leiden(gg, objective_function = "modularity", weights = w, n_iterations = -1)
      list(memb = cm$membership, q = igraph::modularity(gg, cm$membership, weights = w))
    }),
    walktrap = {  # deterministic, no seed needed
      cm <- igraph::cluster_walktrap(gg, weights = w, steps = 4)
      list(membership = igraph::membership(cm),
           Q = igraph::modularity(gg, igraph::membership(cm), weights = w))
    },
    stop("unknown community method: ", method_name)
  )
}

COMMUNITY_METHODS <- c("louvain", "leiden", "walktrap")

# igraph::rewire() keeps each node's degree but DROPS edge attributes (a
# check on the real graphs found every null edge weight was NA, so
# the previous null Q was 0 for every method). So the weights are re-attached
# here: the observed weight multiset is randomly shuffled over the rewired
# edges. Node DEGREES are preserved exactly; node STRENGTHS are not (a hub's
# weights are re-drawn from the whole graph's weight distribution).
rewire_weighted <- function(gg) {
  g_r <- tryCatch(igraph::rewire(gg, with = igraph::keeping_degseq(niter = igraph::gsize(gg) * 10)),
                  error = function(e) NULL)
  if (is.null(g_r)) return(NULL)
  igraph::E(g_r)$weight <- sample(igraph::E(gg)$weight)
  g_r
}

run_modularity_null_test <- function(g, method_name, n_perm = 50) {
  g_use <- igraph::as.undirected(g, mode = "collapse", edge.attr.comb = "sum")
  obs_Q <- run_community_method(g_use, method_name)$Q

  null_Q <- vapply(seq_len(n_perm), function(i) {
    set.seed(10000 + i)  # reproducible rewiring (community runs set their own seeds)
    g_r <- rewire_weighted(g_use)
    if (is.null(g_r)) return(NA_real_)
    res <- tryCatch(run_community_method(g_r, method_name), error = function(e) NULL)
    if (is.null(res)) return(NA_real_)
    res$Q
  }, numeric(1))
  null_Q <- null_Q[!is.na(null_Q)]

  tibble(
    obs_Q = obs_Q,
    null_mean_Q = if (length(null_Q) > 0) mean(null_Q) else NA_real_,
    null_sd_Q   = if (length(null_Q) > 0) stats::sd(null_Q) else NA_real_,
    n_perm_used = length(null_Q),
    p_value = if (length(null_Q) > 0) (sum(null_Q >= obs_Q) + 1) / (length(null_Q) + 1) else NA_real_
  )
}

N_PERM_MODULARITY <- 1000  # 

graphs_1D <- graphs  # already restricted to the minusA19 variants at the top of the script

modularity_null_rows <- list()
for (variant in names(graphs_1D)) {
  for (method_name in COMMUNITY_METHODS) {
    res <- tryCatch(run_modularity_null_test(graphs_1D[[variant]], method_name, n_perm = N_PERM_MODULARITY),
                     error = function(e) {
                       warning("1D: ", variant, "/", method_name, " failed: ", conditionMessage(e))
                       tibble(obs_Q = NA_real_, null_mean_Q = NA_real_, null_sd_Q = NA_real_, n_perm_used = 0L, p_value = NA_real_)
                     })
    modularity_null_rows[[length(modularity_null_rows) + 1]] <- res %>% mutate(variant = variant, method = method_name)
  }
  message("  ", variant, ": modularity null test done.")
}
modularity_null_tbl <- bind_rows(modularity_null_rows) %>% select(variant, method, everything()) %>%
  mutate(method = factor(method, levels = COMMUNITY_METHODS))
openxlsx::write.xlsx(modularity_null_tbl %>% mutate(method = as.character(method)),
                     file.path(OUT_DIR_D, "modularity_null_test.xlsx"), overwrite = TRUE)

p_mod_null <- modularity_null_tbl %>%
  mutate(sig = p_value < 0.05) %>%
  ggplot(aes(x = method, y = obs_Q)) +
  geom_col(aes(fill = sig), color = "black") +
  geom_errorbar(aes(ymin = null_mean_Q - null_sd_Q, ymax = null_mean_Q + null_sd_Q), width = 0.2, color = "grey30") +
  geom_point(aes(y = null_mean_Q), shape = 4, size = 2) +
  facet_wrap(~variant) +
  scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70"), name = "p < 0.05 vs. null") +
  labs(x = NULL, y = "Modularity Q (bar = observed; x = null mean +/- SD)",
       title = "Observed modularity vs. degree/weight-preserving null model (minusA19)") +
  theme_prism_like() + theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR_D, "modularity_null_test.pdf"), p_mod_null, width = 12, height = 9)

message("1D done.")