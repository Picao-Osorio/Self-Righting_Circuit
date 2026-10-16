#!/usr/bin/env Rscript
# ap_lr_connectivity_analysis.R
#
# Anterior-Posterior (axial) and Ipsi/Contra (lateral) connectivity
# polarity analysis -- R port of ap_lr_connectivity_analysis.py, plus the
# "overall" pyramid plots, the "along the body axis" heatmaps, and the
# segment-trend statistical tests, all in one script.
#
# INPUT (both required, flat "detailed" files -- one row per source-neuron /
# partner-neuron edge):
#     ALL_SOURCE_NEURONS_downstream_groups_detailed.xlsx
#     ALL_SOURCE_NEURONS_upstream_groups_detailed.xlsx
#
# USAGE:
#     Rscript ap_lr_connectivity_analysis.R <downstream.xlsx> <upstream.xlsx> <output_dir>
#
# OUTPUT (written to <output_dir>):
#     ap_lr_summary_by_group.csv        -- AP + LR % synapses per group x direction (overall)
#     ap_lr_by_segment.csv              -- net AP / LR polarity per group x direction x segment
#     trend_stats_results.csv           -- ordinal/logistic trend tests along the body axis, FDR-corrected
#     overall_axial_polarity.png        -- mirrored bar chart, upstream left / downstream right
#     overall_lateral_polarity.png      -- mirrored bar chart, upstream left / downstream right
#     axial_polarity_by_segment.png     -- heatmap, downstream + upstream facets, with significance stars
#     lateral_polarity_by_segment.png   -- heatmap, downstream + upstream facets
#
# COLOR CONVENTION (consistent between the overall bars and the body-axis
# heatmaps, and matching ap_lr_connectivity_analysis.py's group definitions):
#     Axial:   Anterior = teal (#1D9E75), Same segment = gray (#888780), Posterior = coral (#D85A30)
#     Lateral: Ipsilateral = pink (#D4537E), Contralateral = purple (#7F77DD)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(MASS)   # polr() for the ordinal (proportional-odds) trend test
  library(svglite) # svg output for ggsave()
})

args <- commandArgs(trailingOnly = TRUE)

# Default file names used when the script is run with NO arguments -- e.g.
# pressing "Source" in RStudio, which does not pass command-line args at all.
# Edit these three lines if your files are named/located differently, or just
# call the script from a terminal with explicit arguments instead (see USAGE
# above): that always overrides these defaults.
DEFAULT_DOWNSTREAM <- "ALL_SOURCE_NEURONS_downstream_groups_detailed_A19fl.xlsx"
DEFAULT_UPSTREAM   <- "ALL_SOURCE_NEURONS_upstream_groups_detailed_A19fl.xlsx"
DEFAULT_OUTPUT_DIR <- "output"

if (length(args) == 0) {
  message("No command-line arguments supplied (this is normal if you ran the ",
          "script via RStudio's Source button rather than from a terminal). ",
          "Falling back to default file names in the current working directory:\n",
          "  downstream : ", DEFAULT_DOWNSTREAM, "\n",
          "  upstream   : ", DEFAULT_UPSTREAM, "\n",
          "  output dir : ", DEFAULT_OUTPUT_DIR, "\n",
          "Current working directory: ", getwd())
  if (!file.exists(DEFAULT_DOWNSTREAM) || !file.exists(DEFAULT_UPSTREAM)) {
    stop("Default files not found in the working directory. Either:\n",
         "  (a) set your working directory to the folder containing the two .xlsx files ",
         "(setwd(...) or RStudio 'Session > Set Working Directory > To Source File Location'), or\n",
         "  (b) run from a terminal with explicit arguments:\n",
         "      Rscript ap_lr_connectivity_analysis_v4.R <downstream.xlsx> <upstream.xlsx> <output_dir>")
  }
  downstream_path <- DEFAULT_DOWNSTREAM
  upstream_path   <- DEFAULT_UPSTREAM
  out_dir         <- DEFAULT_OUTPUT_DIR
} else if (length(args) == 3) {
  downstream_path <- args[1]
  upstream_path   <- args[2]
  out_dir         <- args[3]
} else {
  stop("Usage: Rscript ap_lr_connectivity_analysis_v4.R <downstream.xlsx> <upstream.xlsx> <output_dir>")
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------------
## 1. Canonical source-neuron-group mapping (mirrors 'Matched group(s)')
## ---------------------------------------------------------------------
CANONICAL_GROUP <- c(
  "A01c" = "A01c", "A02o" = "A02o (Wave-1)", "wave" = "A02o (Wave-1)",
  "A05q" = "A05q", "A14a" = "A14a", "A18j" = "A18j",
  "A23g" = "A23g", "A26f" = "A26f", "A27k" = "A27k",
  "basin_2" = "A09a (Basin-2)", "basin_3" = "A09g (Basin-3)", "basin_4" = "A09c (Basin-4)",
  "dnb" = "A09l (DnB)", "goro" = "Goro", "LT1" = "LT1", "LT2" = "LT2", "TePn05" = "TePn05",
  "ddaC" = "mdIV", "v'ada" = "mdIV", "vdaB" = "mdIV",
  "lch5-1" = "chordotonals", "lch5-2/4" = "chordotonals", "lch5-3" = "chordotonals", "lch5-5" = "chordotonals",
  "ddaC?" = "mdIV", "vdaB?" = "mdIV"
)

# Order used throughout the figures (matches the requested neuron ordering)
GROUP_ORDER <- c("mdIV", "chordotonals", "A09a (Basin-2)", "A09g (Basin-3)", "A09c (Basin-4)",
                  "A05q", "A23g", "Goro", "TePn05", "A02o (Wave-1)", "A09l (DnB)", "A27k",
                 "A14a", "A01c", "A18j", "A26f", "LT2", "LT1")

SEGMENT_ORDER <- setNames(c(1, 2, 3, 4:12), c("T1", "T2", "T3", paste0("A", 1:9)))
SEGMENTS <- names(SEGMENT_ORDER)

## ---------------------------------------------------------------------
## 2. Colors -- SAME palette used for the overall bars and the heatmaps
## ---------------------------------------------------------------------
COL_ANTERIOR <- "#1D9E75"   # teal
COL_SAME     <- "#888780"   # gray
COL_POSTERIOR <- "#D85A30"  # coral
COL_IPSI     <- "#D4537E"   # pink
COL_CONTRA   <- "#7F77DD"   # purple

## ---------------------------------------------------------------------
## 3. Load + classify edges
## ---------------------------------------------------------------------
seg_pos <- function(seg) unname(SEGMENT_ORDER[toupper(trimws(seg))])

ap_category <- function(src_seg, partner_seg) {
  si <- seg_pos(src_seg); pi <- seg_pos(partner_seg)
  ifelse(is.na(si) | is.na(pi), "Unresolved",
         ifelse(pi < si, "Anterior", ifelse(pi > si, "Posterior", "Same segment")))
}

lr_category <- function(src_side, partner_side) {
  s <- toupper(trimws(src_side)); p <- toupper(trimws(partner_side))
  ok <- s %in% c("L", "R") & p %in% c("L", "R") & !is.na(s) & !is.na(p)
  ifelse(!ok, "Unresolved", ifelse(s == p, "Ipsilateral", "Contralateral"))
}

load_direction <- function(path, direction) {
  df <- read_excel(path)
  if (direction == "downstream") {
    partner_seg <- df$Downstream_neuron_Segment
    partner_side <- df$Downstream_neuron_Side
    partner_catmaid <- df$Downstream_neuron_CATMAID
    partner_name <- df$Downstream_neuron_name
  } else {
    partner_seg <- df$Upstream_neuron_Segment
    partner_side <- df$Upstream_neuron_Side
    partner_catmaid <- df$Upstream_neuron_CATMAID
    partner_name <- df$Upstream_neuron_name
  }
  out <- tibble(
    source_neuron_CATMAID = df$source_neuron_Neuron_CATMAID,
    source_neuron_name = df$source_neuron_Neuron_name,
    source_segment = df$source_neuron_Segment,
    source_side = df$source_neuron_Side,
    partner_CATMAID = partner_catmaid,
    partner_name = partner_name,
    partner_segment = partner_seg,
    partner_side = partner_side,
    synapses = df$`#Synapses`,
    matched_group = df$`Matched group(s)`,
    direction = direction
  )
  out$AP_category <- ap_category(out$source_segment, out$partner_segment)
  out$LR_category <- lr_category(out$source_side, out$partner_side)
  out
}

down <- load_direction(downstream_path, "downstream")
up   <- load_direction(upstream_path, "upstream")
edges <- bind_rows(down, up)
edges$source_group <- unname(CANONICAL_GROUP[edges$source_neuron_name])
if (any(is.na(edges$source_group))) {
  missing <- unique(edges$source_neuron_name[is.na(edges$source_group)])
  stop(paste("Unmapped source_neuron_name(s), add to CANONICAL_GROUP:", paste(missing, collapse = ", ")))
}
edges$source_group <- factor(edges$source_group, levels = GROUP_ORDER)
edges$source_seg_pos <- seg_pos(edges$source_segment)

cat(sprintf("Loaded %d edges (%d downstream, %d upstream)\n", nrow(edges), nrow(down), nrow(up)))

## ---------------------------------------------------------------------
## 4. Overall summary: % synapses per group x direction x category
## ---------------------------------------------------------------------
pct_summary <- function(data, cat_col, categories) {
  data %>%
    filter(.data[[cat_col]] != "Unresolved") %>%
    group_by(source_group, direction, category = .data[[cat_col]]) %>%
    summarise(synapses = sum(synapses), .groups = "drop") %>%
    group_by(source_group, direction) %>%
    mutate(pct = 100 * synapses / sum(synapses)) %>%
    ungroup() %>%
    complete(source_group, direction, category = categories, fill = list(synapses = 0, pct = NA))
}

ap_overall <- pct_summary(edges, "AP_category", c("Anterior", "Same segment", "Posterior"))
lr_overall <- pct_summary(edges, "LR_category", c("Ipsilateral", "Contralateral"))

write.csv(
  bind_rows(mutate(ap_overall, polarity = "axial"), mutate(lr_overall, polarity = "lateral")),
  file.path(out_dir, "ap_lr_summary_by_group_woutA19fl.csv"), row.names = FALSE
)

## ---------------------------------------------------------------------
## 5. Along the body axis: net polarity per group x direction x segment
## ---------------------------------------------------------------------
net_by_segment <- function(data, cat_col, pos_level, neg_level) {
  data %>%
    filter(.data[[cat_col]] != "Unresolved") %>%
    mutate(segment = factor(source_segment, levels = SEGMENTS)) %>%
    filter(!is.na(segment)) %>%
    group_by(source_group, direction, segment) %>%
    summarise(
      total = sum(synapses),
      pos = sum(synapses[.data[[cat_col]] == pos_level]),
      neg = sum(synapses[.data[[cat_col]] == neg_level]),
      n_edges = n(),
      .groups = "drop"
    ) %>%
    mutate(net = 100 * (pos - neg) / total) %>%
    complete(source_group, direction, segment = factor(SEGMENTS, levels = SEGMENTS))
}

ap_by_seg <- net_by_segment(edges, "AP_category", "Posterior", "Anterior")
lr_by_seg <- edges %>%
  filter(LR_category != "Unresolved") %>%
  mutate(segment = factor(source_segment, levels = SEGMENTS)) %>%
  filter(!is.na(segment)) %>%
  group_by(source_group, direction, segment) %>%
  summarise(total = sum(synapses), ipsi = sum(synapses[LR_category == "Ipsilateral"]),
            n_edges = n(), .groups = "drop") %>%
  mutate(net = 100 * ipsi / total - 50) %>%
  complete(source_group, direction, segment = factor(SEGMENTS, levels = SEGMENTS))

write.csv(bind_rows(mutate(ap_by_seg, polarity = "axial"), mutate(lr_by_seg, polarity = "lateral")),
          file.path(out_dir, "ap_lr_by_segment_woutA19fl.csv"), row.names = FALSE)

## ---------------------------------------------------------------------
## 6. Trend tests along the body axis (edge-level, NOT synapse-weighted --
##    independence assumption: a synapse is not an independent trial, an
##    edge/connection is)
## ---------------------------------------------------------------------
test_ap_trend <- function(sub) {
  d <- sub %>% filter(AP_category != "Unresolved", !is.na(source_seg_pos))
  n_segs <- n_distinct(d$source_seg_pos)
  cat_ok <- sum(table(d$AP_category) >= 2)
  if (nrow(d) < 15 || n_segs < 3 || cat_ok < 2) {
    return(list(n_edges = nrow(d), n_segments = n_segs, coef = NA, p = NA, note = "insufficient data"))
  }
  d$AP_ord <- factor(d$AP_category, levels = c("Anterior", "Same segment", "Posterior"), ordered = TRUE)
  result <- tryCatch({
    fit <- polr(AP_ord ~ source_seg_pos, data = d, Hess = TRUE)
    ct <- coef(summary(fit))
    coef_val <- ct["source_seg_pos", "Value"]
    t_val <- ct["source_seg_pos", "t value"]
    p_val <- 2 * pnorm(-abs(t_val))
    list(n_edges = nrow(d), n_segments = n_segs, coef = coef_val, p = p_val, note = "")
  }, error = function(e) {
    list(n_edges = nrow(d), n_segments = n_segs, coef = NA, p = NA, note = paste("model failed:", conditionMessage(e)))
  })
  result
}

test_lr_trend <- function(sub) {
  d <- sub %>% filter(LR_category != "Unresolved", !is.na(source_seg_pos))
  n_segs <- n_distinct(d$source_seg_pos)
  tab <- table(d$LR_category)
  if (nrow(d) < 15 || n_segs < 3 || length(tab) < 2 || min(tab) < 2) {
    return(list(n_edges = nrow(d), n_segments = n_segs, coef = NA, p = NA, note = "insufficient data"))
  }
  d$y <- as.integer(d$LR_category == "Contralateral")
  result <- tryCatch({
    fit <- glm(y ~ source_seg_pos, data = d, family = binomial)
    s <- summary(fit)$coefficients
    list(n_edges = nrow(d), n_segments = n_segs,
         coef = s["source_seg_pos", "Estimate"], p = s["source_seg_pos", "Pr(>|z|)"], note = "")
  }, error = function(e) {
    list(n_edges = nrow(d), n_segments = n_segs, coef = NA, p = NA, note = paste("model failed:", conditionMessage(e)))
  })
  result
}

combos <- edges %>% distinct(source_group, direction)
trend_rows <- lapply(seq_len(nrow(combos)), function(i) {
  g <- combos$source_group[i]; d <- combos$direction[i]
  sub <- edges %>% filter(source_group == g, direction == d)
  ap_res <- test_ap_trend(sub)
  lr_res <- test_lr_trend(sub)
  tibble(source_group = g, direction = d,
         AP_n_edges = ap_res$n_edges, AP_n_segments = ap_res$n_segments,
         AP_coef = ap_res$coef, AP_p = ap_res$p, AP_note = ap_res$note,
         LR_n_edges = lr_res$n_edges, LR_n_segments = lr_res$n_segments,
         LR_coef = lr_res$coef, LR_p = lr_res$p, LR_note = lr_res$note)
})
trend_df <- bind_rows(trend_rows)
trend_df$AP_q <- p.adjust(trend_df$AP_p, method = "BH")
trend_df$LR_q <- p.adjust(trend_df$LR_p, method = "BH")
write.csv(trend_df, file.path(out_dir, "trend_stats_results_woutA19fl.csv"), row.names = FALSE)

cat("\nSignificant AXIAL trends (FDR q < 0.05):\n")
print(trend_df %>% filter(AP_q < 0.05) %>% dplyr::select(source_group, direction, AP_coef, AP_p, AP_q))
cat("\nSignificant LATERAL trends (FDR q < 0.05):\n")
sig_lr <- trend_df %>% filter(LR_q < 0.05) %>% dplyr::select(source_group, direction, LR_coef, LR_p, LR_q)
if (nrow(sig_lr) == 0) cat("  none\n") else print(sig_lr)

## ---------------------------------------------------------------------
## 7. PLOTS -- all saved as .svg (vector, editable in Illustrator/Inkscape)
## ---------------------------------------------------------------------

## 7a. Overall -- mirrored bar chart, upstream (left, negative) / downstream (right, positive)
plot_pyramid <- function(overall_df, categories, colors, filename, title) {
  df <- overall_df %>%
    mutate(pct = ifelse(is.na(pct), 0, pct),
           category = factor(category, levels = categories),
           signed_pct = ifelse(direction == "upstream", -pct, pct))
  p <- ggplot(df, aes(x = source_group, y = signed_pct, fill = category)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = setNames(colors, categories)) +
    scale_y_continuous(labels = function(v) paste0(abs(v), "%"), limits = c(-101, 100)) +
    scale_x_discrete(limits = rev(GROUP_ORDER)) +
    labs(title = title, x = NULL, y = "% synapses  (upstream \u2190  |  \u2192 downstream)", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")
  ggsave(file.path(out_dir, filename), p, device = "svg", width = 9, height = 8)
  p
}

plot_pyramid(ap_overall, c("Anterior", "Same segment", "Posterior"),
             c(COL_ANTERIOR, COL_SAME, COL_POSTERIOR),
             "overall_axial_polarity.svg", "Axial polarity -- upstream (left) vs downstream (right)")

plot_pyramid(lr_overall, c("Ipsilateral", "Contralateral"),
             c(COL_IPSI, COL_CONTRA),
             "overall_lateral_polarity.svg", "Lateral polarity -- upstream (left) vs downstream (right)")

## 7b. Along the body axis -- heatmap, faceted by direction, with classical
##     tiered significance stars (* q<0.05, ** q<0.01, *** q<0.001, FDR-corrected).
##     show_sig = TRUE  -> adds the "sig." annotation column (stats version)
##     show_sig = FALSE -> plain heatmap, no stats column (clean version)
sig_stars <- function(q) {
  dplyr::case_when(
    is.na(q) ~ "",
    q < 0.001 ~ "***",
    q < 0.01  ~ "**",
    q < 0.05  ~ "*",
    TRUE ~ ""
  )
}

plot_heatmap <- function(by_seg_df, trend_df, q_col, low_col, high_col, filename, title, legend_lab,
                          show_sig = TRUE) {
  x_levels <- if (show_sig) c(SEGMENTS, "sig") else SEGMENTS
  df <- by_seg_df %>% mutate(segment = factor(as.character(segment), levels = x_levels))

  if (show_sig) {
    # one extra "sig" column per group x direction showing the classical stars
    sig_df <- trend_df %>%
      mutate(stars = sig_stars(.data[[q_col]])) %>%
      dplyr::select(source_group, direction, stars) %>%
      mutate(segment = factor("sig", levels = x_levels), net = NA_real_)
    plot_df <- bind_rows(df %>% mutate(stars = NA_character_), sig_df)
  } else {
    plot_df <- df %>% mutate(stars = NA_character_)
  }

  p <- ggplot(plot_df, aes(x = segment, y = source_group, fill = net)) +
    geom_tile(data = filter(plot_df, segment != "sig"), color = "white", linewidth = 0.5) +
    geom_text(data = filter(plot_df, segment != "sig"),
              aes(label = ifelse(is.na(net), "", round(net))), size = 2.6, color = "black") +
    facet_wrap(~direction, ncol = 2) +
    scale_fill_gradient2(low = low_col, mid = COL_SAME, high = high_col, midpoint = 0,
                          limits = c(-100, 100), na.value = "white", name = legend_lab) +
    scale_y_discrete(limits = rev(GROUP_ORDER)) +
    labs(title = title, x = "source segment", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 0))

  if (show_sig) {
    p <- p +
      geom_text(data = filter(plot_df, segment == "sig"),
                aes(label = stars), size = 3.4, fontface = "bold", color = "black") +
      scale_x_discrete(labels = c(setNames(SEGMENTS, SEGMENTS), sig = "sig.")) +
      labs(caption = "sig. column: trend test along the body axis, FDR-corrected -- * q<0.05  ** q<0.01  *** q<0.001")
  }

  ggsave(file.path(out_dir, filename), p, device = "svg",
         width = if (show_sig) 11.5 else 10.5, height = 8)
  p
}

# stats-annotated versions
plot_heatmap(ap_by_seg, trend_df, "AP_q", COL_ANTERIOR, COL_POSTERIOR,
             "axial_polarity_by_segment.svg",
             "Axial polarity along the body axis",
             "net (post-ant)", show_sig = TRUE)

plot_heatmap(lr_by_seg, trend_df, "LR_q", COL_CONTRA, COL_IPSI,
             "lateral_polarity_by_segment.svg",
             "Lateral polarity along the body axis",
             "net (ipsi-50)", show_sig = TRUE)

# plain versions, no stats column
plot_heatmap(ap_by_seg, trend_df, "AP_q", COL_ANTERIOR, COL_POSTERIOR,
             "axial_polarity_by_segment_nostats.svg",
             "Axial polarity along the body axis",
             "net (post-ant)", show_sig = FALSE)

plot_heatmap(lr_by_seg, trend_df, "LR_q", COL_CONTRA, COL_IPSI,
             "lateral_polarity_by_segment_nostats.svg",
             "Lateral polarity along the body axis",
             "net (ipsi-50)", show_sig = FALSE)

cat(sprintf("\nDone. Outputs written to: %s\n", normalizePath(out_dir)))
