# Data acquisition for the fitting chapter.
#
# Targets:
#   make            same as `make fetch`
#   make fetch      run R script to write *.tsv from outbreaks pkg
#   make clean      remove fetched TSVs
#
# The BMJ 1978 PDF (under source/) is committed to the repo — URLs of
# 50-year-old single-page mirrors rot, and the citations against it are
# load-bearing. Don't fetch it dynamically.
#
# Requirements:
#   R (with the `outbreaks` package — auto-installed on first run)
#
# Outputs:
#   influenza_1978_school.tsv      wide canonical: time, date,
#                                  in_bed, convalescent
#   in_bed.tsv                     narrow: time, in_bed
#                                  (consumed by the fit pipeline)

TSVS := influenza_1978_school.tsv in_bed.tsv

.PHONY: all fetch clean

all: fetch

fetch: $(TSVS)

$(TSVS): fetch.R
	Rscript fetch.R

clean:
	rm -f $(TSVS)
