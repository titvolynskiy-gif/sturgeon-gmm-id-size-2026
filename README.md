# Sturgeon Cranial Bone Shape Analysis

Geometric morphometric analysis of cranial bones in four native sturgeon
species (*Huso huso*, *Acipenser gueldenstaedtii*, *A. ruthenus*,
*A. stellatus*), *A. baerii*, and two hybrids (Bester, Sevbel).

Code accompanying the article:
> Tyt Volynskyi1,2, Svitozar Davydenko2, Oleksandr Kovalchuk1,3,4,*

1 National Museum of Natural History, National Academy of Sciences of Ukraine, 15 Bohdana Khmelnytskoho St, Kyiv 01054, Ukraine

2 Schmalhausen Institute of Zoology, National Academy of Sciences of Ukraine, 15 Bohdana Khmelnytskoho St, Kyiv 01054, Ukraine

3 Department of Palaeozoology, Faculty of Biological Sciences, University of Wrocław, 21 Sienkiewicza St, Wrocław 50-335, Poland

4 Department of Biology and Biology Teaching Methodology, Faculty of Natural Sciences and Geography, A. S. Makarenko Sumy State Pedagogical University, 87 Romenska St, Sumy 40002, Ukraine
2026. *Title of the article*. Journal, DOI: [add DOI]

## What the script does

1. Generalized Procrustes Analysis (GPA) for the six 2D-digitised bones;
   reads the five bones that were already registered in 3D.
2. Averages left/right sides per individual and runs a per-bone PCA.
3. Tests taxonomic shape signal (`procD.lm`), for all four species and
   for *Acipenser* only.
4. Tests allometry (shape ~ log(centroid size) × species).
5. Tests bilateral (directional/fluctuating) asymmetry for paired bones.
6. Runs pairwise PERMANOVA between species for every bone.
7. Builds Ward.D2 dendrograms from per-bone PC scores (4-species core set
   and the full 7-taxon set).
8. Produces all figures used in the article (Figs. 1–4 + per-bone PCA
   scatterplots).

## Repository structure

```
.
├── data/                              # input landmark tables (see below)
│   ├── all_landmarks_2D.csv
│   └── all_landmarks_3D.csv
├── scripts/
│   └── morphometrics_analysis.R
├── LICENSE
├── CITATION.cff
├── .gitignore
└── README.md
```

Running the script creates an `output/` folder (`tables/` + `figures/`)
next to the data. That folder is **not** tracked in git (see
`.gitignore`) — it is fully reproducible from the code + data.

## Requirements

- R ≥ 4.2
- packages: `geomorph`, `vegan`, `ggplot2`, `ggdendro`

```r
install.packages(c("geomorph", "vegan", "ggplot2", "ggdendro"))
```

## How to run

1. Clone the repository.
2. Open `scripts/morphometrics_analysis.R` and set `DATA_DIR` to the
   folder with the two CSV files (e.g. `"data"` if you run the script
   from the repository root).
3. Run the script top to bottom.

## Input data format

`all_landmarks_2D.csv` / `all_landmarks_3D.csv` — one row per digitised
specimen/side, long format:

| column | description |
|---|---|
| `bone` | bone name |
| `species` | species / hybrid name |
| `fish` | individual fish ID |
| `side` | `L`/`R` for paired bones, blank otherwise |
| `n_lm` | number of landmarks for that bone |
| `LM1_X, LM1_Y, (LM1_Z)...` | landmark coordinates |
| `centroid` | centroid size (3D file only; for 2D it comes from GPA) |

## Output

- `output/tables/` — CSV results (sample sizes, PCA coverage, taxonomic
  signal, allometry, bilateral symmetry, pairwise PERMANOVA, dendrogram
  objects, `sessionInfo.txt`).
- `output/figures/` — all article figures + per-bone PCA scatterplots.

## License

Code: MIT License (see `LICENSE`).
Data (`data/*.csv`): [state your chosen data license here, e.g. CC-BY 4.0].

## Citation

If you use this code, please cite the article above. You can additionally
cite the archived code release itself (see `CITATION.cff`; Zenodo DOI:
[add after archiving the repository on Zenodo]).
