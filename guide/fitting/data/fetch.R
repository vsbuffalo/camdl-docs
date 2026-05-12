#!/usr/bin/env Rscript
# Fetch the 1978 English boarding-school influenza dataset from the
# R `outbreaks` package and emit two TSV artefacts.
#
# Source dataset:
#   `outbreaks::influenza_england_1978_school` — daily counts of boys
#   confined to bed and convalescent at an English boarding school
#   during a 14-day H1N1 outbreak in January 1978.
#
# Original publication:
#   Anonymous (1978). "Influenza in a boarding school." British Medical
#   Journal 1(6112): 587. PMID: 626866.
#   PDF mirrored at data/source/bmj_1978_influenza_boarding_school.pdf
#
# Outputs (relative to the data/ directory this script lives in):
#   influenza_1978_school.tsv  — wide canonical export: time, date,
#                                in_bed, convalescent.
#   in_bed.tsv                 — narrow file consumed by the fit
#                                pipeline (time, in_bed only).

# --- Resolve script location so paths work from any CWD ---
this_script <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=--file=).+", args, perl = TRUE))
  if (length(m) > 0) return(normalizePath(m[1]))
  if (sys.nframe() > 0) return(normalizePath(sys.frame(1)$ofile))
  NA_character_
}
DATA_DIR <- dirname(this_script())
if (is.na(DATA_DIR)) DATA_DIR <- "."

# --- Install outbreaks on demand ---
if (!requireNamespace("outbreaks", quietly = TRUE)) {
  message("installing R package: outbreaks")
  install.packages("outbreaks", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages(library(outbreaks))
d <- influenza_england_1978_school

# Day 0 = first observation (1978-01-22)
d$time <- as.integer(d$date - d$date[1])

wide <- data.frame(
  time         = d$time,
  date         = format(d$date, "%Y-%m-%d"),
  in_bed       = d$in_bed,
  convalescent = d$convalescent
)

# Wide canonical
wide_path <- file.path(DATA_DIR, "influenza_1978_school.tsv")
write.table(wide, file = wide_path, sep = "\t",
            row.names = FALSE, quote = FALSE)

# Narrow for the fit pipeline (camdl reads `time` + `in_bed`)
narrow_path <- file.path(DATA_DIR, "in_bed.tsv")
write.table(wide[, c("time", "in_bed")], file = narrow_path,
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("wrote %s (%d rows)\n",
            wide_path, nrow(wide)))
cat(sprintf("wrote %s (%d rows)\n",
            narrow_path, nrow(wide)))
cat(sprintf("source: outbreaks::influenza_england_1978_school   "))
cat(sprintf("(R outbreaks pkg version %s)\n",
            packageVersion("outbreaks")))
