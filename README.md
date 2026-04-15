# camdl-book

Documentation, tutorials, and vignettes for [camdl](https://github.com/idm-tools/camdl).
Builds the site at https://vincebuffalo.org/camdl/docs/.

## Structure

```
guide/          Getting started — language tutorial, first fit
language/       Language spec, data model (synced from camdl repo)
inference/      Conceptual guide, fitting pipeline, debugging
vignettes/      Worked examples — spatial inference, sampler comparison, etc.
reference/      Feature reference, specs (synced from camdl repo), design philosophy
models/         .camdl files used by tutorials and vignettes
data/           Datasets (synthetic and real)
configs/        Fit configs (.toml) for tutorials and vignettes
figures/        Hand-made figures (tikz sources, PNGs)
```

## Building

Requires [Quarto](https://quarto.org/) and a working `camdl` binary on PATH.
The camdl repo must be checked out as a sibling directory (or set `CAMDL_REPO`).

```bash
make render      # sync specs from camdl repo, then quarto render
make preview     # live preview
make deploy      # rsync to server
```

## Synced specs

Some reference chapters are pulled from `camdl/docs/` at build time via
`make sync`. The canonical copy of these files lives in the
[camdl repo](https://github.com/idm-tools/camdl) — do not edit them here.

| Book page | Source |
|-----------|--------|
| `language/spec.qmd` | `camdl/docs/camdl-language-spec.md` |
| `language/data-spec.qmd` | `camdl/docs/camdl-data-spec.md` |
| `reference/inference-spec.qmd` | `camdl/docs/camdl-inference-spec.md` |
| `reference/ir-spec.qmd` | `camdl/docs/compartmental-ir-spec.md` |
| `reference/runtimes.qmd` | `camdl/docs/runtimes.md` |
| `inference/debugging.qmd` | `camdl/docs/debugging.md` |

## The rule

**Spec or developer reference → [`camdl/docs/`](https://github.com/idm-tools/camdl/tree/main/docs).**
Would a code change require updating this doc in the same commit? It stays there.

**Tutorial, guide, vignette, or figure → here.**
Is this about how to *think about* or *use* camdl? It lives in this repo.
