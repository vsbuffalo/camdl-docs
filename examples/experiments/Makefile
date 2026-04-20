.PHONY: data clean

DATA_DIR = data

data: $(DATA_DIR)/boarding_school_flu.tsv

$(DATA_DIR)/boarding_school_flu.tsv: scripts/fetch_boarding_school.R
	@mkdir -p $(DATA_DIR)
	Rscript $<

clean:
	rm -rf $(DATA_DIR)
