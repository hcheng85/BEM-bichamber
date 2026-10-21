# Hydrodynamic R50/gap50 Fokker--Planck validation

This directory contains only the standalone analysis completed for the current
round. Raw simulation MAT files, processed MAT files, and extracted NPZ caches
are deliberately excluded.

## Geometry convention

The archive calls the available ensemble `R50/gap50`. Its stored physical
parameters give `gap/R=1.0` and `gap/(2R)=0.5`. No `R50/gap25` hydrodynamic
ensemble, which would literally give `gap/R=0.5`, was available.

## Reproduction

Create an environment using `requirements.txt`, then extract compact features:

```bash
python scripts/extract_gap50_hydro_features.py \
  --raw-dir /path/to/BEM_data \
  --processed-dir /path/to/Processed_data \
  --output-dir /path/to/feature_cache
```

Run the validation without modifying the source data:

```bash
python scripts/analyze_gap50_hydro_fp.py \
  --processed-dir /path/to/Processed_data \
  --feature-dir /path/to/feature_cache \
  --result-dir /path/to/regenerated_results
```

The scripts expect the 20 raw `...R50_num100_wid50_n*_BEM.mat` files and their
processed `event_histogram_...R50_gap50_n*_BEM.mat` counterparts. See
`configuration.json` for the exact split, parameters, lags, thresholds, and
smoothing lengths. `REPORT.md` gives the scientific conclusions and
`summary.csv` provides the principal numerical validation metrics.
