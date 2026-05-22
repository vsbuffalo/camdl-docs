#!/usr/bin/env Rscript
# Hagelloch 1861 measles outbreak -> daily incidence TSV for camdl.
#
# Source: outbreaks::measles_hagelloch_1861 (Pfeilsticker 1863; Oesterle 1992).
# 188 children, one German village, a single school-class introduction.
# We aggregate prodrome-onset dates to a daily case-incidence series and anchor
# t = 0 at 1861-10-01 (safely before the index introduction), so the seed time
# tau is a positive parameter the model can estimate.
suppressMessages(library(outbreaks))

d <- measles_hagelloch_1861
T0 <- as.Date("1861-10-01")

onset <- as.Date(d$date_of_prodrome)
day <- as.integer(onset - T0)
counts <- tabulate(day + 1L, nbins = max(day) + 1L)   # day 0..max
inc <- data.frame(time = 0:(length(counts) - 1L), cases = counts)

out <- "hagelloch_daily.tsv"   # run from guide/fitting/seed-timing/data/
write.table(inc, out, sep = "\t", row.names = FALSE, quote = FALSE)

# --- provenance / validation anchor -------------------------------------------
# The epidemiological index is the earliest-onset case. The introduction (its
# infection) precedes the prodrome by the measles incubation period (~10-11 d),
# so the seed time tau should land in mid-October — about a week and a half
# before the first observed onset. (The SIR has no explicit latent stage, so
# tau is the effective infection-onset time, not a literal calendar instant.)
idx <- which.min(day)
cat(sprintf("wrote %s: %d days, %d total cases (t=0 = %s)\n",
            out, nrow(inc), sum(counts), T0))
cat(sprintf("index case ID %d: prodrome %s (day %d), rash %s\n",
            d$case_ID[idx], as.character(onset[idx]), min(day),
            as.character(as.Date(d$date_of_rash[idx]))))
cat(sprintf("epidemic peak: day %d (%d cases); slow early phase days %d-%d then takeoff\n",
            which.max(counts) - 1L, max(counts), min(day), min(day) + 21L))
