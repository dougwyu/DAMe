#!/usr/bin/env Rscript
# Turn a per-PCR, sequence-level vsearch table into an occupancy-model survey table.
#
# Usage: Rscript perpcr_to_occupancy.R table.tsv clusters.uc PCRinfo.txt survey.tsv
#
#   table.tsv    vsearch --otutabout output: one row per DAMe-passed sequence,
#                one column per PCR (from convert --per-pcr labels)
#   clusters.uc  vsearch --cluster_* --uc output for the same sequences
#   PCRinfo.txt  from dame convert --per-pcr
#   survey.tsv   written: one row per PCR (PCRinfo columns), one column per OTU
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)

# vsearch drops ";size=N;" style annotations from table row names but keeps them
# in .uc files, so compare labels without annotations
strip_annotations <- function(label) sub(";.*$", "", label)

seq_table <- read.delim(args[1], check.names = FALSE) |>
  rename(sequence = 1) |>
  mutate(sequence = strip_annotations(sequence))

# .uc: S rows are centroids (V9 = label), H rows are members (V9 = label, V10 = centroid)
uc <- read.delim(args[2], header = FALSE, stringsAsFactors = FALSE)
membership <- bind_rows(
  uc |> filter(V1 == "S") |> transmute(sequence = V9, centroid = V9),
  uc |> filter(V1 == "H") |> transmute(sequence = V9, centroid = V10)
) |>
  transmute(sequence = strip_annotations(sequence),
            otu = paste0("OTU_", strip_annotations(centroid)))

unassigned <- setdiff(seq_table$sequence, membership$sequence)
if (length(unassigned) > 0) {
  stop("sequences in ", args[1], " missing from ", args[2], ": ",
       paste(head(unassigned, 5), collapse = ", "))
}

pcr_info <- read.delim(args[3], stringsAsFactors = FALSE)

# Sum sequence rows within OTUs
otu_table <- seq_table |>
  inner_join(membership, by = "sequence") |>
  group_by(otu) |>
  summarise(across(-sequence, sum), .groups = "drop")

# One row per PCR, one column per OTU
otu_by_pcr <- otu_table |>
  pivot_longer(-otu, names_to = "pcr_id", values_to = "count") |>
  pivot_wider(names_from = otu, values_from = count)
otu_cols <- setdiff(names(otu_by_pcr), "pcr_id")

extra <- setdiff(otu_by_pcr$pcr_id, pcr_info$pcr_id)
if (length(extra) > 0) {
  stop("PCRs in ", args[1], " missing from ", args[3], ": ",
       paste(head(extra, 5), collapse = ", "))
}

# Keep every PCR in PCRinfo; PCRs with no column in table.tsv had no matched reads
survey <- pcr_info |>
  left_join(otu_by_pcr, by = "pcr_id") |>
  mutate(across(all_of(otu_cols), \(x) replace_na(x, 0L)))

write.table(survey, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
