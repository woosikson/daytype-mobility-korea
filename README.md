# Day-type-resolved mobility matrices by age — code

Code for the Data Descriptor *"Beyond monthly commuting flows: mobility matrices by age for weekdays, weekends and holidays in South Korea"* (Son, Wi, Nah & Jung). It builds, for each month from January 2018 to December 2023 and each day type (weekday, weekend, holiday), origin–destination matrices for 250 districts × 15 age groups divided into commute, non-commute and non-move.

- Data: https://doi.org/10.5281/zenodo.23092377

## Method in brief

1. Initial guess from public data: census commuting shares × registered population with a radiation model (commute), KTDB purpose-specific OD and Personal Travel Survey participation (non-commute).
2. Weekend and holiday initial guesses from Time Use Survey coefficients and θ = (φ, ψ, χ).
3. Iterative proportional fitting (IPF) to the registered population (origins) and the 15:00 mobile-phone de facto population averaged by day type (destinations); decomposition into commute, non-commute and non-move.
4. Bayesian calibration of θ in September 2022 by grid integration; the posterior median is used for all 72 months.

See the Methods of the paper for equations and details.

## Layout

| Path | Content |
|---|---|
| `run.jl`, `src/*.jl` | All computation (Julia). Writes `data/*.csv` and `data/od/*.csv.gz`. |
| `plot.R`, `_common.R` | Figures (R, ggplot2). Reads only `data/*.csv`; writes `figs/*.pdf`, `*.png`, `*.eps`. |
| `flowchart.tex` | Figure 1 (TikZ). |
| `Project.toml`, `Manifest.toml` | Julia environment. |
| `input/` | Input data (not in this repository; see below). |
| `contents/` | Time Use Survey report PDF (not in this repository; see below). |

## Inputs

1. Download `inputs_public.zip` from the Zenodo data record and unzip it into `input/`.
2. The following inputs cannot be redistributed; obtain them from the providers and place them in `input/` (file names as listed; formats are described in the data record README):
   - `defacto_daytype_2018-2023.csv` — monthly day-type averages of the daily mobile-phone de facto population (KT Corporation);
   - `od_corrected_202209.csv` — mobile-phone day/night residence OD, September 2022 (KT Corporation; validation only);
   - `nc_od_agegroup.csv` — KTDB OD by purpose (shopping, leisure, other), 2019 (https://www.ktdb.go.kr);
   - `nc_person_row_factor.csv`, `person_factors.csv` — Personal Travel Survey 2021 participation (https://www.ktdb.go.kr);
   - `visitor_202209.csv`, `visitor_sido_202209.csv` — Korea Tourism Data Lab visitor data (validation only);
   - `sgg_polygons.csv` — district boundary polygons (SGIS; maps only).
3. Place the 2024 Time Use Survey report `제1-1권 생활시간량편(표1-5).pdf` (Ministry of Data and Statistics) in `contents/`. The code extracts Tables 1-2 and 2-2 with `pdftotext` (poppler), which must be installed.

Without the restricted inputs the code cannot be run end to end; the released matrices in the data record are its output.

## Run

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia -t 12 --project=. run.jl          # about 26 minutes on 12 threads (Apple M2 Max); deterministic, no random numbers
Rscript plot.R                           # about 3 minutes
pdflatex flowchart.tex                   # Figure 1
```

Software: Julia 1.12.4 (IPF via `ProportionalFitting.jl`; versions fixed in `Manifest.toml`); R 4.5 with ggplot2, dplyr, tidyr, readr, scales, cowplot, viridisLite, ragg, gtable, latex2exp and (optional, for evenly spaced text in vector output) showtext.

## License

The code in this repository is released under the [MIT License](LICENSE). The data records (the 174 OD matrices and accompanying files on Zenodo) are released separately under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) license. Restricted inputs (KT, KTDB and Korea Tourism Data Lab data) are not covered by either license and remain subject to the terms of their providers.

## Citation

If you use this code or the data, please cite the Data Descriptor (Beyond monthly commuting flows: mobility matrices by age for weekdays, weekends and holidays in South Korea) and the data record (DOI: 10.5281/zenodo.23092377).

