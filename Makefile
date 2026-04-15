CAMDL_REPO ?= ../camdl
REMOTE     := vincebuffalo:apps/vincebuffalo/camdl/docs
SITE       := _site

SPECS := \
	language/spec.qmd:$(CAMDL_REPO)/docs/camdl-language-spec.md:Language\ Specification \
	language/data-spec.qmd:$(CAMDL_REPO)/docs/camdl-data-spec.md:Data\ Model \
	reference/inference-spec.qmd:$(CAMDL_REPO)/docs/camdl-inference-spec.md:Inference\ Specification \
	reference/ir-spec.qmd:$(CAMDL_REPO)/docs/compartmental-ir-spec.md:IR\ Specification \
	reference/runtimes.qmd:$(CAMDL_REPO)/docs/runtimes.md:Simulation\ Backends \
	inference/debugging.qmd:$(CAMDL_REPO)/docs/debugging.md:Debugging

.PHONY: sync render preview deploy clean

sync:
	@for entry in $(SPECS); do \
		dest=$${entry%%:*}; rest=$${entry#*:}; \
		src=$${rest%%:*}; title=$${rest#*:}; \
		printf -- '---\ntitle: "%s"\n---\n\n' "$$title" > $$dest; \
		echo '<!-- DO NOT EDIT — synced from camdl/docs/ via make sync -->' >> $$dest; \
		echo '' >> $$dest; \
		cat "$$src" >> $$dest; \
		echo "  $$src -> $$dest"; \
	done

render: sync
	@command -v camdl >/dev/null 2>&1 || { echo "error: camdl not on PATH (install from ../camdl via make install)"; exit 1; }
	@command -v camdlc >/dev/null 2>&1 || { echo "error: camdlc not on PATH (install from ../camdl via make install)"; exit 1; }
	uv run quarto render

preview:
	uv run quarto preview

deploy: render
	rsync -avz --delete $(SITE)/ $(REMOTE)

clean:
	rm -rf $(SITE) .quarto
