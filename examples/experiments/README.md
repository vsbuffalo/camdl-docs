# Simulation and Experiments — example files

Example models, parameters, and data for the
[Simulation and Experiments](../../guide/experiments.qmd) chapter.

## Data

Run `make data` to fetch datasets used in this chapter:

```bash
make data
```

This downloads the 1978 English boarding school influenza dataset from the
R `outbreaks` package and writes it to `data/boarding_school_flu.tsv`.

### Boarding school influenza (1978)

- **Source:** Anonymous (1978). "Influenza in a boarding school."
  *British Medical Journal*, 1:578.
- **Population:** 763 boys at a boarding school in North England.
- **Outbreak:** H1N1 influenza, 22 January – 4 February 1978.
  512 of 763 boys were confined to bed (67% attack rate).
- **Data:** Daily prevalence — number of boys confined to bed and
  number convalescent, over 14 days.
- **R package:** `outbreaks::influenza_england_1978_school`

## Models

- `sirv.camdl` — SIR with vaccination (V compartment), used for
  intervention and scenario comparison examples.
- `boarding_school.camdl` — SIR model for the 1978 boarding school
  outbreak, with informed priors on influenza parameters.

## Running

```bash
# Simulate baseline and vaccination scenarios
camdl simulate sirv.camdl --scenario baseline --seed 42
camdl simulate sirv.camdl --scenario with_vaccination --seed 42

# Prior predictive check
camdl simulate boarding_school.camdl --draws prior -n 200 --scenario baseline --seed 42
```
