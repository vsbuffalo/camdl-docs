#!/usr/bin/env Rscript
# Fetch the 1978 English boarding school influenza dataset
# from the R 'outbreaks' package and write to TSV.
#
# Source: Anonymous (1978) "Influenza in a boarding school."
#         British Medical Journal, 1:578.
#         763 boys, 512 infected, H1N1 influenza.
#
# The dataset is daily prevalence (number confined to bed)
# over 14 days (1978-01-22 to 1978-02-04).

if (!requireNamespace("outbreaks", quietly = TRUE)) {
  install.packages("outbreaks", repos = "https://cloud.r-project.org")
}

library(outbreaks)
d <- influenza_england_1978_school

# Add a day column (0-indexed from outbreak start)
d$day <- as.integer(d$date - d$date[1])

# Write TSV (output path relative to script location)
script_dir <- getwd()  # Makefile runs from examples/experiments/
out <- file.path(script_dir, "data", "boarding_school_flu.tsv")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
write.table(d[, c("day", "date", "in_bed", "convalescent")],
            file = out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote", out, "\n")
