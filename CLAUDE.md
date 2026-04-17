# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

R package `opencanopy` that adapts the **Open-Canopy** pre-trained models (UNet/SMP and PVTv2/timm, originally trained on SPOT 6-7 imagery at 1.5 m) to run on **IGN BD ORTHO®** aerial photography at 0.20 m. The end-to-end goal is a canopy height model (CHM) raster for an arbitrary AOI.

End-user workflow is typically: provide an AOI as a GeoPackage → download IGN ortho via WMS → resample → run inference in Python (via reticulate) → export CHM TIFFs.

The package language (code, identifiers, comments, user-facing messages) is **French**. Preserve that convention.

## Common commands

All commands assume the repo root as working directory.

```bash
# End-to-end smoke test (CLI entry point)
Rscript inst/scripts/test_pipeline_aoi.R                            # uses data/aoi.gpkg, pvtv2
Rscript inst/scripts/test_pipeline_aoi.R path/to/aoi.gpkg unet      # custom AOI + model

# Regenerate NAMESPACE + man/ from roxygen comments
Rscript -e 'devtools::document()'

# Install the package locally
R CMD INSTALL .
# or: Rscript -e 'devtools::install()'

# R CMD check (data/ and outputs/ are in .Rbuildignore)
R CMD build . && R CMD check opencanopy_*.tar.gz

# Tests (testthat is declared in Suggests but no tests/ directory is checked in yet)
Rscript -e 'devtools::test()'
```

### Running from an R session (interactive)

The four `R/` files are designed to be `source()`'d directly — they are not just package code, they are also callable as scripts. See `README.md` for the high-level API (`pipeline_aoi_to_chm()`, `predict_chm_from_ign()`, `compute_ndvi()`, …).

### Python/conda environment

Inference requires a conda env named **`open_canopy`** with: `torch`, `torchvision`, `numpy`, `rasterio`, `segmentation-models-pytorch`, `timm`. `reticulate` auto-detects it via `conda_list()`, `CONDA_PREFIX`, or a pre-activated Python on `PATH` (see `inst/scripts/test_pipeline_aoi.R` for the detection chain).

## Architecture

### Four-stage R pipeline (`R/0{1..4}_*.R`)

The numbering is the data-flow order, not just filename sorting:

1. **`download_open_canopy.R`** — acquisition. Hugging Face dataset download (native `hfhub` preferred, Python `huggingface_hub` fallback), IGN WMS download (`data.geopf.fr/wms-r`), loading `.jp2` ortho, resampling ortho IGN (0.20 m) → SPOT grid (1.5 m) with `aggregate(fact = 8)`.
2. **`analyse_open_canopy.R`** — spectral indices (NDVI, GNDVI, SAVI, NDWI) from the IRC ortho (bands PIR/R/V), vegetation masks, multi-panel plots, zonal stats, NDVI × CHM cross-tabulation.
3. **`prediction_open_canopy.R`** — inference. Sets up the conda env, downloads pretrained checkpoints from HF (auto-discovers the checkpoint filename via the HF API in `find_checkpoint_name()`), tiles/stitches for inference, calls Python via `reticulate`, upsamples the 1.5 m CHM back to 0.20 m.
4. **`pipeline_aoi_to_chm.R`** — orchestration. `pipeline_aoi_to_chm("aoi.gpkg")` glues the above: AOI → Lambert-93 reprojection → tiled WMS download (`WMS_MAX_PX = 4096`) → RVB+IRC combination into 4-band (R, G, B, NIR) → aggregation → Python inference → mosaic → CHM TIFFs + indices (NDVI/GNDVI/SAVI/NDWI) + masked CHM + recap PDF (3×3 panels).

`NAMESPACE` is grouped by source file (see comments) — keep that grouping when adding exports.

### R ↔ Python bridge

- `reticulate::use_condaenv("open_canopy", required = TRUE)` is the canonical way; the test script has fallback detection for RStudio launched outside the conda env.
- UNet checkpoints load directly via `segmentation_models_pytorch`.
- PVTv2 requires extra code from the upstream repo. This is handled in two ways:
  1. Preferred: **embedded standalone module** at `inst/python/timmnet_standalone.py` — vendored minimal subset of `fajwel/Open-Canopy` (`timmNet`, `SimpleSegmentationHead`, `infer_output`, `set_first_layer`). This avoids cloning the upstream repo.
  2. Fallback: `download_open_canopy_src()` clones the upstream repo on demand (see `OPEN_CANOPY_SRC` config in `04_*.R`).

### Model contract — do not break these invariants

These are stated at the top of `R/pipeline_aoi_to_chm.R` and the recent commit history (e.g. `9687201 "remove incorrect input normalization"`) confirms they have already bitten:

- **Input:** 4 channels in order **R, G, B, NIR** (SPOT 6-7 order).
- **No normalization.** The model expects **raw pixel values**, not 0–1 or mean/std normalized inputs. Earlier versions normalized and produced wrong predictions.
- **Targets stored as decimetres / 10** in the Open-Canopy dataset; the model outputs metres directly.
- **Checkpoint layout:** PyTorch Lightning, weights prefixed with `net.seg_model.` — strip this prefix when loading into a bare `smp.Unet` / `timmNet`.
- **Resampling factor is 8** (0.20 m → 1.5 m ≈ ×7.5, rounded up to 8). Pre-agg with `terra::aggregate(fact = 8, fun = "mean")`; post-disagg with `disagg(fact = 8)` if you want CHM at native IGN resolution.

### Geospatial conventions

- All AOI inputs are reprojected to **Lambert-93 (EPSG:2154)** before any WMS call — IGN layers are served in L93.
- WMS layer names are millésime-aware: `ign_layer_name()` in `04_*.R` appends the year (`ORTHOIMAGERY.ORTHOPHOTOS2024`, `ORTHOIMAGERY.ORTHOPHOTOS.IRC.2024`); `NULL` means "most recent".
- Any AOI > ~25 km² will take a long time (multi-tile WMS + tiled inference); the test script warns but does not abort.

### Outputs

`outputs/` is gitignored. Expect: `ortho_rvb.tif`, `ortho_irc.tif`, `ndvi.tif`, `chm_predicted_1_5m.tif`, `chm_predicted_0_2m.tif`, `resultats_aoi.pdf`.

## Notes for edits

- The code base is flat (no sub-modules in `R/`). Function-level docs are roxygen2; regenerating with `devtools::document()` rewrites both `NAMESPACE` and `man/`.
- HF auth: `Sys.setenv(HF_TOKEN = "hf_...")` is read by both the `hfhub` R path and the Python fallback.
- Recent work has focused on (a) switching HF downloads to the R-native `hfhub` path (commit `343847c`), (b) Windows symlink warnings in `.resolve_hf_path()` (`fc6427d`), and (c) downsampling before RStudio preview (`d4634b1`) — keep these fixes in mind when touching the download / display code.
