# opencanopy 0.1.1

Maintenance release consolidating a session of downstream-integration
fixes discovered while wiring the pipeline into the `nemetonshiny`
Shiny app.

### New Features

* **Phase + tile progress callback** — `pipeline_aoi_to_chm()`,
  `download_ortho_for_aoi()`, `download_ign_tiled()` and
  `run_inference()` now accept an optional `progress_callback` that
  fires with a structured list at each pipeline boundary:
    * `list(type = "phase", step = <load_aoi | download_ortho |
      setup_python | download_model | inference | export>, total = 5L,
      model = <opt>)`
    * `list(type = "tile_phase_start", prefix, n_tiles)` +
      `list(type = "tile", prefix, idx, n_tiles, bbox)` per WMS tile
    * `list(type = "inference_phase_start", n_tiles, model)` +
      `list(type = "inference_tile", idx, n_tiles, model, tile_name)`
  Used by nemetonshiny to paint a live status line instead of an
  8-minute silent wait. Backward compatible: `NULL` callback is a
  no-op.

### Bug Fixes

* **Robust Python discovery** — `setup_python()` no longer trusts the
  first `py_config()$python` it finds. It now walks a prioritised
  list (`RETICULATE_PYTHON` → currently-bound Python if it matches →
  canonical miniforge3 / miniconda3 / anaconda3 / mambaforge install
  paths → `CONDA_PREFIX` → `conda_list()` as last resort) and
  aborts with an explicit error when reticulate is already bound
  to a Python that is *not* the `open_canopy` env, telling the user
  exactly which `Sys.setenv(RETICULATE_PYTHON = …)` to run.
* **`do.call(merge, predictions)`** — now qualified as
  `do.call(terra::merge, unname(predictions))` in `run_inference()`
  and `download_ign_tiled()`. On AOIs large enough to need > 1
  inference tile, the named list fed to the bare `merge()` used to
  dispatch to `base::merge.default()` which errored with
  "l'argument x est manquant, avec aucune valeur par défaut" after
  5+ minutes of successful inference.
* **`fs::file_delete()` qualified** everywhere — the unqualified
  `file_delete()` call in the cleanup path of `download_ign_tiled()`
  hit `NAMESPACE` gaps (only `fs::dir_ls` / `dir_create` were
  imported) and failed with "impossible de trouver la fonction
  file_delete", leaving 56 orphan tile rasters behind.
* **Tile-sidecar cleanup** — the post-mosaic cleanup now also
  removes the `*.tif.aux.xml` and `*.tif.ovr` GDAL sidecars written
  by terra alongside each tile, so project caches no longer drift
  towards 100+ orphan `.xml` files.

### Documentation

* **`setup_python()` docstring and error messages** clarified: the
  R `hfhub` package replaces only the Python `huggingface_hub`
  sub-module for the checkpoint download step. Python (torch +
  rasterio + segmentation-models-pytorch + timm) is always required
  for inference. The previous wording hinted that `install.packages
  ('hfhub')` could bypass Python entirely, which was misleading.

# opencanopy 0.1.0

Premier release du package `opencanopy`.

## Point d'entree principal

- `pipeline_aoi_to_chm(aoi_path)` : pipeline complet AOI (GeoPackage) -> ortho IGN -> CHM predit + indices spectraux + PDF recapitulatif.

## Fonctionnalites

- Telechargement automatique des orthophotos IGN (BD ORTHO RVB + IRC, 0.20 m) via le WMS Geoplateforme, avec decoupage par tuiles pour les AOI > 4096 px.
- Telechargement des checkpoints Open-Canopy (`unet`, `pvtv2`) depuis Hugging Face via `hfhub` (R natif) ou Python `huggingface_hub` en fallback.
- Inference PyTorch via `reticulate` sur un environnement conda `open_canopy` (UNet/SMP ou PVTv2/timm reconstruit via module Python embarque `inst/python/timmnet_standalone.py`).
- Calcul de 4 indices spectraux depuis l'IRC : NDVI, GNDVI, SAVI, NDWI.
- CHM masque vegetation (`chm_vegetation_0_2m.tif`) via seuils NDVI/NDWI configurables.
- Reechantillonnage automatique 0.20 m -> 1.5 m pour l'inference, puis disaggregation bilineaire 1.5 m -> 0.20 m sur la sortie.
- Visualisation recapitulative (`resultats_aoi.pdf`) en grille 3x3 : orthos, indices, masque, CHM brut et masque.
- Support des millesimes IGN (`millesime_ortho`, `millesime_irc`).

## Parametres exposes dans `pipeline_aoi_to_chm()`

- `model_name` : `"pvtv2"` (defaut, plus precis) ou `"unet"` (plus rapide).
- `ndvi_threshold` (defaut 0.25) et `ndwi_threshold` (defaut 0.20) pour le masque vegetation.
- `millesime_ortho`, `millesime_irc` pour cibler un millesime specifique.

## Dependances

- Packages R : `terra`, `sf`, `httr2`, `jsonlite`, `curl`, `fs`, `reticulate`, `hfhub` (suggere).
- Environnement Python : `torch`, `torchvision`, `numpy`, `rasterio`, `segmentation-models-pytorch`, `timm`, `torchmetrics`, `omegaconf`, `hydra-core`, `lightning`, `einops`.
