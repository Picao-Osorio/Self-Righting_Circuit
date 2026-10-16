## ---------------------------------------------------------------------
## Shared data-prep for the FULLY segment/side-resolved connectivity
## analysis -- both the SOURCE neuron instance and the DOWNSTREAM/TARGET
## neuron instance are resolved to their own (group, segment, side).
##
## the previous version pooled every downstream target to its group-level sum (using the Totals
## sheet's pre-summed "Sum <target> downstream #Synapses" columns). What
## is actually wanted is the per-edge detail: e.g. A01c_a1l has 28
## synapses onto LT1_a1r specifically, out of A01c_a1l's Total #Synapses
## of 193 -> proportion 28/193 = 0.145. That edge-level detail lives in
## the "Filtered_downstream" sheet, not "Totals" -- so "Totals" is now
## used ONLY as the source of each source-instance's Total #Synapses
## (the proportion denominator); every actual connection (numerator) is
## read from Filtered_downstream.
##
## Sourced by both analyze_connectivity_matrix_segment_side.R and
## gephi_style_network_segment_side.R.
##
## Produces THREE resolution variants, each as a pair of matrices (raw #Synapses / arcsin(sqrt(
## proportion))):
##   1. side x side  -- mat_side_side_raw / mat_side_side_arcsin
##        rows = source instance   (Group_SegSide)
##        cols = target instance   (Group_SegSide)
##        Fully resolved on both axes -- this is the "A01c_a1l -> LT1_a1r"
##        level of detail.
##   2. seg x seg    -- mat_seg_seg_raw / mat_seg_seg_arcsin
##        rows = source group+segment (L/R averaged)
##        cols = target group+segment (L/R averaged)
##        Both sides collapsed, mirroring the original group-level
##        L/R-averaging on both axes.
##   3. seg x side   -- mat_seg_side_raw / mat_seg_side_arcsin
##        rows = source group+segment (L/R averaged, as in variant 2)
##        cols = target instance      (Group_SegSide, as in variant 1)
##        Only the axis that was collapsed before (source) is collapsed;
##        the newly-resolved target axis is left alone.
##
## Averaging convention (unchanged from the group-level pipeline): raw
## #Synapses is averaged with a plain mean; the arcsin transform is
## applied to the MEAN proportion, never averaged after transforming
## (mean/median and arcsin(sqrt(.)) do not commute) -- see
## /projects/.../conventions-and-decisions.md.
##
## ---- Handling of malformed/missing-position rows -----------------------
## Originally 62 of 1213 Filtered_downstream rows had a missing/invalid
## downstream segment or side. Manually curation of the workbook
## (ALL_SOURCE_NEURONS_downstream_groups_summary-mod_flagged-mod.xlsx):
## 
##   (a) 46 rows, downstream = "TePn05 classIV-related projection C": the
##       segment was recoverable and filled in as A8 (side was already
##       present) -- these now resolve normally like any other row.
##   (b) 16 rows remain genuinely unresolvable -- 5 "vdaB?" (Matched
##       group(s) = mdIV, segment+side both blank) and 11 corrupted
##       lch5-2/4 / lch5-5 entries (Matched group(s) = chordotonals,
##       segment and/or side missing or invalid).these 16 rows are EXCLUDED from the
##       analysis entirely (not bucketed into an "unresolved" column/node
##       as the previous version did) -- their #Synapses do not count
##       toward any target instance, though the source instance's
##       Total_Synapses denominator (from Totals) is unaffected.
## -----------------------------------------------------------------------

required_pkgs <- c("readxl", "dplyr", "tidyr", "tibble", "stringr")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
invisible(lapply(required_pkgs, library, character.only = TRUE))

DATA_FILE <- "ALL_SOURCE_NEURONS_downstream_groups_summary-mod_flagged-mod.xlsx"

# ---- 0. Fixed lookup tables (shared with the group-level v3/v9 scripts,
# and with the earlier group-pooled version of this segment/side script) --

segment_levels <- c(paste0("T", 1:3), paste0("A", 1:9))

hierarchy_rank <- c(
  "mdIV" = 1, "chordotonals" = 1,
  "A09a (Basin-2)" = 2, "A09g (Basin-3)" = 2, "A09c (Basin-4)" = 2,
  "A05q" = 3, "A23g" = 3, "Goro" = 4, "TePn05" = 4, "A02o (Wave-1)" = 4, "A19f" = 5,
  "A19l" = 6, "A14a" = 6, "A01c" = 6, "A18j" = 6, "A26f" = 6, "A09l (DnB)" = 4, "A27k" = 5,
  "LT1" = 7, "LT2" = 7
)

neuron_color_map <- c(
  "LT1"             = "#FF9300",
  "LT2"             = "#0000FF",
  "A26f"            = "#996633",
  "A18j"            = "#996633",
  "A01c"            = "#996633",
  "A14a"            = "#996633",
  "A19l"            = "#996633",
  "A09l (DnB)"      = "#FF2600",
  "A27k"            = "#FF66CC",
  "A19f"            = "#FF66CC",
  "A02o (Wave-1)"   = "#800000",
  "TePn05"          = "#400040",
  "Goro"            = "#2C3E50",
  "A23g"            = "#FFD479",
  "A05q"            = "#FFD479",
  "A09c (Basin-4)"  = "#008F00",
  "A09g (Basin-3)"  = "#008F00",
  "A09a (Basin-2)"  = "#008F00",
  "chordotonals"    = "#004080",
  "mdIV"            = "#004080"
)

neuron_type_map <- c(
  "mdIV" = "sensory", "chordotonals" = "sensory",
  "A09a (Basin-2)" = "interneuron", "A09g (Basin-3)" = "interneuron", "A09c (Basin-4)" = "interneuron",
  "A05q" = "interneuron", "A23g" = "interneuron",
  "Goro" = "interneuron", "TePn05" = "interneuron",
  "A02o (Wave-1)" = "interneuron", "A09l (DnB)" = "interneuron",
  "A27k" = "interneuron", "A19f" = "interneuron", "A26f" = "interneuron",
  "A18j" = "interneuron", "A01c" = "interneuron", "A14a" = "interneuron",
  "A19l" = "interneuron",
  "LT1" = "motoneuron", "LT2" = "motoneuron"
)

get_colors_by_group <- function(groups) {
  cols <- neuron_color_map[groups]
  if (any(is.na(cols))) {
    message("No color defined for: ", paste(unique(groups[is.na(cols)]), collapse = ", "),
            " -- defaulting to grey60.")
    cols[is.na(cols)] <- "grey60"
  }
  unname(cols)
}

get_circuit_rank <- function(groups) {
  r <- hierarchy_rank[groups]
  if (any(is.na(r))) {
    message("No circuit-hierarchy rank defined for: ", paste(unique(groups[is.na(r)]), collapse = ", "),
            " -- defaulting to mid-rank (3).")
    r[is.na(r)] <- 3
  }
  unname(r)
}

get_type <- function(groups) {
  ty <- neuron_type_map[groups]
  if (any(is.na(ty))) {
    message("No neuron type defined for: ", paste(unique(groups[is.na(ty)]), collapse = ", "),
            " -- defaulting to interneuron (hexagon).")
    ty[is.na(ty)] <- "interneuron"
  }
  unname(ty)
}

# Segment index. Anything outside T1-T3/A1-A9 falls through to "sorted
# last" with a message -- a safety net (no such value is expected to reach
# here any more: rows with an unresolvable segment/side are excluded
# upstream in fd_parsed, not passed through).
seg_index <- function(segs) {
  idx <- match(segs, segment_levels)
  if (any(is.na(idx))) {
    not_in_list <- unique(segs[is.na(idx)])
    not_in_list <- not_in_list[not_in_list != "unresolved"]
    if (length(not_in_list) > 0) {
      message("Segment(s) not in the expected T1-T3/A1-A9 list, sorted last: ",
              paste(not_in_list, collapse = ", "))
    }
    idx[is.na(idx)] <- length(segment_levels) + 1
  }
  idx
}

# Side index: L before R. "none" (an axis that has no side dimension at
# all, e.g. a segment-only L/R-averaged node in the network plots) sorts
# before both; "unresolved" is kept only as a safety net (see seg_index
# above -- no row/node should carry it any more).
side_index_map <- c(none = 0, L = 1, R = 2, unresolved = 3)

# ---- 1. Totals sheet: ONLY used for each source instance's Total #Synapses
# (the proportion denominator). Parsing mirrors the earlier group-pooled
# script, but nothing else from Totals (its 20 "Sum ... downstream
# #Synapses" columns) is used any more -- those were the group-pooled sums
# this rebuild moves away from. ---------------------------------------

totals_raw <- read_excel(DATA_FILE, sheet = "Totals")
names(totals_raw) <- str_trim(names(totals_raw))
totals_raw <- totals_raw %>% mutate(`Source neuron` = str_trim(`Source neuron`))

if (!is.numeric(totals_raw[["Total #Synapses"]])) {
  bad <- totals_raw[["Total #Synapses"]][is.na(suppressWarnings(as.numeric(totals_raw[["Total #Synapses"]]))) &
                                            !is.na(totals_raw[["Total #Synapses"]])]
  if (length(bad) > 0) {
    message("Totals$`Total #Synapses` has non-numeric entries: ", paste(unique(bad), collapse = ", "),
            " (", length(bad), " cell(s)) -- coercing to NA then treating as 0.")
  }
}
totals_raw <- totals_raw %>%
  mutate(`Total #Synapses` = suppressWarnings(as.numeric(`Total #Synapses`))) %>%
  mutate(`Total #Synapses` = replace_na(`Total #Synapses`, 0))

# Raw source-name prefix -> canonical group name, AS USED IN THE TOTALS
# SHEET (which already pools ddaC/v'ada/vdaB into "mdIV_..." rows, and
# chordotonal components into "chordotonal(s)_..." rows). The Totals
# sheet's own source-row prefix is actually PLURAL ("chordotonals_a1l",
# etc. -- confirmed from the "could not be parsed/mapped" warning, which
# was firing because this map only had the singular "chordotonal" key).
# Both spellings are kept here (singular, in case an older/mixed workbook
# uses it; plural, which is what the current Totals sheet actually
# contains) so a re-export either way still resolves to the same
# canonical name.  "chordotonal" and "chordotonals" are merged into the
#  single canonical name "chordotonals" (plural) everywhere in this pipeline.

raw_prefix_map_totals <- c(
  "A01c" = "A01c", "A05q" = "A05q", "A14a" = "A14a", "A18j" = "A18j", "A19f" = "A19f",
  "A19l" = "A19l", "A23g" = "A23g", "A26f" = "A26f", "A27k" = "A27k", "TePn05" = "TePn05",
  "basin_2" = "A09a (Basin-2)", "basin_3" = "A09g (Basin-3)", "basin_4" = "A09c (Basin-4)",
  "chordotonal" = "chordotonals", "chordotonals" = "chordotonals", "dnb" = "A09l (DnB)",
  "goro" = "Goro", "mdIV" = "mdIV", "wave" = "A02o (Wave-1)",
  "ddaC?" = "mdIV", "vdaB?" = "mdIV"
)

totals_parsed <- totals_raw %>%
  extract(`Source neuron`, into = c("raw_prefix", "seg_letter", "seg_num", "side"),
          regex = "^(.*)_([AaTt])(\\d+)([LlRr])$", remove = FALSE) %>%
  mutate(Seg = paste0(str_to_upper(seg_letter), seg_num),
         Side = str_to_upper(side))

unmatched_tot <- totals_parsed %>% filter(is.na(raw_prefix) | !(raw_prefix %in% names(raw_prefix_map_totals)))
if (nrow(unmatched_tot) > 0) {
  message("WARNING: ", nrow(unmatched_tot), " Totals row(s) could not be parsed/mapped and are EXCLUDED ",
          "from the Total_Synapses lookup: ", paste(unmatched_tot$`Source neuron`, collapse = ", "))
}

totals_parsed <- totals_parsed %>%
  filter(!is.na(raw_prefix), raw_prefix %in% names(raw_prefix_map_totals)) %>%
  mutate(Group = unname(raw_prefix_map_totals[raw_prefix]))

canonical_totals <- totals_parsed %>%
  group_by(Group, Seg, Side) %>%
  summarise(Total_Synapses = sum(`Total #Synapses`, na.rm = TRUE), .groups = "drop") %>%
  mutate(Row_label_side = paste0(Group, "_", Seg, Side))

total_lookup <- setNames(canonical_totals$Total_Synapses, canonical_totals$Row_label_side)

# ---- 2. Filtered_downstream sheet: the actual per-edge detail ---------
# One row per (source neuron instance, downstream neuron, segment/side).
# This is the sheet that actually carries "A01c_a1l has 28 synapses onto
# LT1_a1r".

fd_raw <- read_excel(DATA_FILE, sheet = "Filtered_downstream")
names(fd_raw) <- str_trim(names(fd_raw))
fd_raw <- fd_raw %>%
  mutate(
    source_neuron_Neuron_CATMAID = str_trim(source_neuron_Neuron_CATMAID),
    `Matched group(s)`           = str_trim(`Matched group(s)`),
    Downstream_neuron_Segment    = str_trim(Downstream_neuron_Segment),
    Downstream_neuron_Side       = str_trim(Downstream_neuron_Side)
  )

if (!is.numeric(fd_raw[["#Synapses"]])) {
  bad <- fd_raw[["#Synapses"]][is.na(suppressWarnings(as.numeric(fd_raw[["#Synapses"]]))) &
                                  !is.na(fd_raw[["#Synapses"]])]
  if (length(bad) > 0) {
    message("Filtered_downstream$`#Synapses` has non-numeric entries: ", paste(unique(bad), collapse = ", "),
            " (", length(bad), " cell(s)) -- coercing to NA then treating as 0.")
  }
}
fd_raw <- fd_raw %>%
  mutate(`#Synapses` = suppressWarnings(as.numeric(`#Synapses`))) %>%
  mutate(`#Synapses` = replace_na(`#Synapses`, 0))

# Source-name prefix -> canonical group, AS USED IN Filtered_downstream
# (which, unlike Totals, lists each sensory COMPONENT as its own row --
# ddaC, ddaC?, v'ada, vdaB, vdaB? separately for mdIV; lch5-1, lch5-2/4,
# lch5-3, lch5-5 separately for chordotonals). ddaC?/vdaB? are folded into
# mdIV per the same 2026-09-23 confirmation used for the Totals sheet.
filtered_prefix_map <- c(
  "A01c" = "A01c", "A05q" = "A05q", "A14a" = "A14a", "A18j" = "A18j", "A19f" = "A19f",
  "A19l" = "A19l", "A23g" = "A23g", "A26f" = "A26f", "A27k" = "A27k", "TePn05" = "TePn05",
  "basin_2" = "A09a (Basin-2)", "basin_3" = "A09g (Basin-3)", "basin_4" = "A09c (Basin-4)",
  "dnb" = "A09l (DnB)", "goro" = "Goro", "wave" = "A02o (Wave-1)",
  "ddaC" = "mdIV", "ddaC?" = "mdIV", "vdaB" = "mdIV", "vdaB?" = "mdIV", "v'ada" = "mdIV",
  "lch5-1" = "chordotonals", "lch5-2/4" = "chordotonals", "lch5-3" = "chordotonals", "lch5-5" = "chordotonals"
)

fd_parsed <- fd_raw %>%
  extract(source_neuron_Neuron_CATMAID, into = c("src_prefix", "seg_letter", "seg_num", "side"),
          regex = "^(.*)_([AaTt])(\\d+)([LlRr])$", remove = FALSE) %>%
  mutate(Seg_src = paste0(str_to_upper(seg_letter), seg_num),
         Side_src = str_to_upper(side))

unmatched_fd <- fd_parsed %>% filter(is.na(src_prefix) | !(src_prefix %in% names(filtered_prefix_map)))
if (nrow(unmatched_fd) > 0) {
  message("WARNING: ", nrow(unmatched_fd), " Filtered_downstream row(s) have a source instance that could not ",
          "be parsed/mapped and are EXCLUDED: ", paste(unmatched_fd$source_neuron_Neuron_CATMAID, collapse = ", "))
}

missing_group_tgt <- fd_parsed %>% filter(is.na(`Matched group(s)`) | `Matched group(s)` == "")
if (nrow(missing_group_tgt) > 0) {
  message("WARNING: ", nrow(missing_group_tgt), " Filtered_downstream row(s) have no `Matched group(s)` value ",
          "(no clean downstream-group identity) and are EXCLUDED.")
}

fd_parsed <- fd_parsed %>%
  filter(!is.na(src_prefix), src_prefix %in% names(filtered_prefix_map),
         !is.na(`Matched group(s)`), `Matched group(s)` != "") %>%
  mutate(Group_src = unname(filtered_prefix_map[src_prefix]),
         Group_tgt = `Matched group(s)`)

# Validate each row's downstream position: BOTH segment and side must be
# present AND valid (segment in T1-T3/A1-A9, side in L/R). Any row that
# still fails this after the manual curation is EXCLUDED from the
# analysis (see header note) -- as of the curated workbook, this should
# only be the 5 "vdaB?" (mdIV) and 11 corrupted lch5-2/4 / lch5-5
# (chordotonals) rows; the 46 TePn05 rows were fixed (segment filled in as
# A8) and now pass this check like any other row.
fd_parsed <- fd_parsed %>%
  mutate(
    seg_tgt_valid  = !is.na(Downstream_neuron_Segment) & Downstream_neuron_Segment %in% segment_levels,
    side_tgt_valid = !is.na(Downstream_neuron_Side) & Downstream_neuron_Side %in% c("L", "R"),
    Position_resolved = seg_tgt_valid & side_tgt_valid
  )

n_excluded <- sum(!fd_parsed$Position_resolved)
if (n_excluded > 0) {
  message("EXCLUDING ", n_excluded, " Filtered_downstream row(s) (",
          sum(fd_parsed$`#Synapses`[!fd_parsed$Position_resolved]),
          " synapses) with a still-missing/invalid downstream segment or side",
          "Affected Matched group(s): ",
          paste(sort(unique(fd_parsed$Group_tgt[!fd_parsed$Position_resolved])), collapse = ", "),
          ". (Source rows: ", paste(fd_parsed$source_neuron_Neuron_CATMAID[!fd_parsed$Position_resolved], collapse = ", "), ")")
}

fd_parsed <- fd_parsed %>%
  filter(Position_resolved) %>%
  mutate(
    Seg_tgt         = Downstream_neuron_Segment,
    Side_tgt        = Downstream_neuron_Side,
    Row_label_side  = paste0(Group_src, "_", Seg_src, Side_src),
    Row_label_seg   = paste0(Group_src, "_", Seg_src),
    Col_label_side  = paste0(Group_tgt, "_", Seg_tgt, Side_tgt),
    Col_label_seg   = paste0(Group_tgt, "_", Seg_tgt)
  )

# ---- 3. Collapse duplicate raw rows sharing the same (source instance,
# target instance) -- e.g. more than one sensory component (ddaC AND
# v'ada) can independently synapse onto the same downstream instance at
# the same segment/side; these are summed, not averaged (same
# same-anatomical-location convention as the group-level source pooling).

edges_side <- fd_parsed %>%
  group_by(Group_src, Seg_src, Side_src, Row_label_side, Row_label_seg,
           Group_tgt, Col_label_side) %>%
  summarise(Synapses = sum(`#Synapses`, na.rm = TRUE), n_pooled = dplyr::n(), .groups = "drop")

n_pooled_combos <- sum(edges_side$n_pooled > 1)
if (n_pooled_combos > 0) {
  message(n_pooled_combos, " (source instance, target instance) edge(s) were pooled (summed) from more than ",
          "one Filtered_downstream row (e.g. two sensory components both synapsing on the same downstream ",
          "instance).")
}

edges_side <- edges_side %>% mutate(Total_Synapses = unname(total_lookup[Row_label_side]))

missing_denom <- edges_side %>% filter(is.na(Total_Synapses))
if (nrow(missing_denom) > 0) {
  message("WARNING: ", nrow(missing_denom), " edge(s) reference a source instance with no Total_Synapses in the ",
          "Totals sheet and are EXCLUDED: ", paste(unique(missing_denom$Row_label_side), collapse = ", "))
  edges_side <- edges_side %>% filter(!is.na(Total_Synapses))
}

edges_side <- edges_side %>%
  mutate(Proportion = ifelse(Total_Synapses > 0, Synapses / Total_Synapses, NA_real_),
         Proportion_arcsin = asin(sqrt(Proportion)))

# Note (not a data error): Totals-sheet source instances that have a
# nonzero Total_Synapses but NO rows at all in Filtered_downstream. 
# 
# Filtered_downstream is deliberately filtered down to the neurons composing the circuit of interest -- so an
# instance like A09g (Basin-3)_A3L (Total #Synapses = 36 in Totals, 0 rows
# in Filtered_downstream) simply has all 36 of its synapses onto
# downstream partners OUTSIDE that circuit, not a spreadsheet gap. Such an
# instance correctly has no row/node in the resolved matrices/networks
# below (there is nothing within the circuit for it to connect to) -- this
# is flagged for transparency only, not as something to fix.
totals_without_edges <- canonical_totals %>%
  filter(Total_Synapses > 0, !(Row_label_side %in% edges_side$Row_label_side))
if (nrow(totals_without_edges) > 0) {
  message(nrow(totals_without_edges), " source instance(s) have a nonzero Total #Synapses in the Totals sheet but ",
          "no rows in Filtered_downstream -- all their synapses go to downstream partners outside the filtered ",
          "circuit of interest, so they have no row/node in the resolved matrices/networks: ",
          paste0(totals_without_edges$Row_label_side, " (Total=", totals_without_edges$Total_Synapses, ")",
                 collapse = ", "))
}

# ---- 4. Row/column metadata (Group/Seg/Side per label, for coloring and
# for the dendrogram-rotation / node-position ranks). Built once here and
# reused by all three variants below. ------------------------------------

row_meta_side <- edges_side %>%
  distinct(Row_label_side, Group = Group_src, Seg = Seg_src, Side = Side_src)

row_meta_seg <- edges_side %>%
  distinct(Row_label_seg, Group = Group_src, Seg = Seg_src)

col_meta_side <- fd_parsed %>%
  distinct(Col_label_side, Group = Group_tgt, Seg = Seg_tgt, Side = Side_tgt)

col_meta_seg <- fd_parsed %>%
  distinct(Col_label_seg, Group = Group_tgt, Seg = Seg_tgt)

# ---- 5. Variant 1: side x side (fully resolved both axes) -------------

mat_side_side_raw <- edges_side %>%
  select(Row_label_side, Col_label_side, Synapses) %>%
  pivot_wider(names_from = Col_label_side, values_from = Synapses, values_fill = 0) %>%
  column_to_rownames("Row_label_side") %>%
  as.matrix()

mat_side_side_arcsin <- edges_side %>%
  select(Row_label_side, Col_label_side, Proportion_arcsin) %>%
  pivot_wider(names_from = Col_label_side, values_from = Proportion_arcsin, values_fill = 0) %>%
  column_to_rownames("Row_label_side") %>%
  as.matrix()

# ---- 6. Variant 2: seg x seg (both sides averaged) --------------------
# mean of raw Synapses across whichever (source-side, target-side) combos
# exist for a given (source Group+Seg, target Group+Seg); mean of the
# FRACTION Proportion likewise, with arcsin(sqrt(.)) applied to that mean
# afterward (never averaged post-transform).

long_seg_seg <- edges_side %>%
  left_join(fd_parsed %>% distinct(Col_label_side, Col_label_seg), by = "Col_label_side") %>%
  group_by(Row_label_seg, Col_label_seg) %>%
  summarise(
    Synapses_mean   = mean(Synapses, na.rm = TRUE),
    Proportion_mean = mean(Proportion, na.rm = TRUE),
    n_combos        = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(Proportion_mean_arcsin = asin(sqrt(Proportion_mean)))

mat_seg_seg_raw <- long_seg_seg %>%
  select(Row_label_seg, Col_label_seg, Synapses_mean) %>%
  pivot_wider(names_from = Col_label_seg, values_from = Synapses_mean, values_fill = 0) %>%
  column_to_rownames("Row_label_seg") %>%
  as.matrix()

mat_seg_seg_arcsin <- long_seg_seg %>%
  select(Row_label_seg, Col_label_seg, Proportion_mean_arcsin) %>%
  pivot_wider(names_from = Col_label_seg, values_from = Proportion_mean_arcsin, values_fill = 0) %>%
  column_to_rownames("Row_label_seg") %>%
  as.matrix()

# ---- 7. Variant 3: seg (source averaged) x side (target resolved) -----

long_seg_side <- edges_side %>%
  group_by(Row_label_seg, Col_label_side) %>%
  summarise(
    Synapses_mean   = mean(Synapses, na.rm = TRUE),
    Proportion_mean = mean(Proportion, na.rm = TRUE),
    n_combos        = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(Proportion_mean_arcsin = asin(sqrt(Proportion_mean)))

mat_seg_side_raw <- long_seg_side %>%
  select(Row_label_seg, Col_label_side, Synapses_mean) %>%
  pivot_wider(names_from = Col_label_side, values_from = Synapses_mean, values_fill = 0) %>%
  column_to_rownames("Row_label_seg") %>%
  as.matrix()

mat_seg_side_arcsin <- long_seg_side %>%
  select(Row_label_seg, Col_label_side, Proportion_mean_arcsin) %>%
  pivot_wider(names_from = Col_label_side, values_from = Proportion_mean_arcsin, values_fill = 0) %>%
  column_to_rownames("Row_label_seg") %>%
  as.matrix()

# ---- 8. Dendrogram-rotation / node-position ranks ----------------------
# Every axis (row or column, side-resolved or seg-only) uses the same
# nesting: axial segment (coarsest) -> side, where the axis has one
# (L before R; "unresolved" trails after both) -> circuit-hierarchy
# (sensory -> ... -> motor, finest). "unresolved" segments/sides always
# sort last via seg_index()/side_index_map's fallback behavior.

rank_side <- function(meta) {
  segi <- seg_index(meta$Seg)
  sidi <- unname(side_index_map[meta$Side])
  circ <- get_circuit_rank(meta$Group)
  segi * 1000 + sidi * 100 + circ
}

rank_seg <- function(meta) {
  segi <- seg_index(meta$Seg)
  circ <- get_circuit_rank(meta$Group)
  segi * 10 + circ
}

rank_lookup_row_side <- setNames(rank_side(row_meta_side), row_meta_side$Row_label_side)
rank_lookup_col_side <- setNames(rank_side(col_meta_side), col_meta_side$Col_label_side)
rank_lookup_row_seg  <- setNames(rank_seg(row_meta_seg),  row_meta_seg$Row_label_seg)
rank_lookup_col_seg  <- setNames(rank_seg(col_meta_seg),  col_meta_seg$Col_label_seg)

reorder_rows_by_rank <- function(mat, rank_lookup) {
  hc <- hclust(dist(mat, method = "euclidean"), method = "average")
  wts <- rank_lookup[hc$labels]
  as.hclust(reorder(as.dendrogram(hc), wts, agglo.FUN = mean))
}

reorder_cols_by_rank <- function(mat, rank_lookup) {
  hc <- hclust(dist(t(mat), method = "euclidean"), method = "average")
  wts <- rank_lookup[hc$labels]
  as.hclust(reorder(as.dendrogram(hc), wts, agglo.FUN = mean))
}

# Group lookup vectors (label -> Group), for coloring row/column labels
# and network nodes by neuron identity.
row_group_side <- setNames(row_meta_side$Group, row_meta_side$Row_label_side)
row_group_seg  <- setNames(row_meta_seg$Group,  row_meta_seg$Row_label_seg)
col_group_side <- setNames(col_meta_side$Group, col_meta_side$Col_label_side)
col_group_seg  <- setNames(col_meta_seg$Group,  col_meta_seg$Col_label_seg)

message("segment_side_data_prep.R: built variant 1 (side x side) ", nrow(mat_side_side_raw), " x ",
        ncol(mat_side_side_raw), ", variant 2 (seg x seg) ", nrow(mat_seg_seg_raw), " x ", ncol(mat_seg_seg_raw),
        ", variant 3 (seg x side) ", nrow(mat_seg_side_raw), " x ", ncol(mat_seg_side_raw), ". Total synapses ",
        "represented: ", sum(edges_side$Synapses), " (Filtered_downstream total: ", sum(fd_raw$`#Synapses`), ").")
