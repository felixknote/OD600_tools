# OD600 Tools

Python scripts for OD600-based growth analysis of CRISPRi knockdown strains in 96-well plates.

**Environment:** `mAIcrobe` conda environment

---

## Scripts

### `01_plate_heatmap_edge_effect.py` — Plate Heatmap & Edge-Effect Analysis

Visualises endpoint OD600 measurements as 96-well plate heatmaps and quantifies the well-position (edge) effect on growth.

**Input:** Plate CSV files (`P<bio>_<tech>.csv`, semicolon-delimited, comma decimal) in the script directory.

**Output:** `plate_analysis/`
- `P1_1.png`, `P1_2.png`, … — individual plate heatmaps
- `average_plate.png` — mean across all plates
- `correlation_od_vs_distance.png` — violin plot of OD600 by edge distance with Pearson r and pairwise Mann-Whitney U tests
- `zone_map.png` — 96-well layout coloured by edge distance (0 = outer edge, 3 = center)

---

### `02_growth_curve.py` — Kinetic Growth Curve Analysis (CRISPRi)

Analyses 20 h kinetic OD600 measurements from 28 CRISPRi target genes (3 sgRNA guides each) across 3 biological replicates.

**Input:**
- 3 kinetic OD600 CSV files from the plate reader
- Plate map CSV (8 x 12, semicolon-delimited)

**Analysis:**
- Per-curve 4PL sigmoidal fits for outlier detection (EC50 +/- 2 SD threshold)
- Mean +/- SD growth curves per sgRNA guide (outlier-excluded)
- Single 4PL fit on the retained mean with EC50 +/- SD annotation
- Overview grid (7x4, 16:9) of all 28 mutants vs WT reference

**Output:** `<data_dir>/growth_curve_analysis/`
```
per_mutant/
  raw/        -- individual replicate curves per gene (30 PNGs)
  summary/    -- mean +/- SD ribbons + 4PL fit per gene (30 PNGs)
overview_grid.png
growth_parameters.csv
growth_parameters_sorted_EC50.xlsx
```

**Growth parameters columns:**
`gene, pathway, bottom, top, ec50_h, ec50_std, hill_slope, max_growth_rate, max_OD600_mean, max_OD600_std, n_curves`

---

## Dependencies

```
numpy  pandas  matplotlib  scipy  openpyxl  tqdm
```

```bash
conda activate mAIcrobe
```
