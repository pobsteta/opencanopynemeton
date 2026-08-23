# opencanopy 1.0.0

Première version majeure. Elle solde le balayage `values()` ouvert par la 0.1.3
et rend rattachable le CRS de tous les rasters produits par le pipeline.

### Bug Fixes

* **Balayage `values()` — les huit sites restants** — la 0.1.3 corrigeait
  `pct_veg` mais laissait en suspens le `grep -n "values(" R/` qu'elle
  réclamait elle-même. Deux sites portaient la même classe de risque :
  `compute_irc_stats()`, appelé sur une **ortho IRC entière** (boucle sur les
  dalles), et `evaluate_predictions()`, qui tenait **trois** couches pleine
  résolution en mémoire simultanément — `prediction`, `reference` et leur
  différence — pour quatre scalaires. Les six autres :
  `normalize_for_model()` (le plus systématique du paquet : chaque tuile de
  chaque prédiction, alors que seul le maximum est utilisé),
  `compute_chm_stats()`, `plot_chm_histogram()`, `plot_canopy_change()` et
  `plot_prediction_comparison()` (deux occurrences). Mesuré à `memfrac = 0.3` :
  `compute_chm_stats()` sur 1,44e8 cellules passe de **4,7 Go à 1,16 Go**
  (85 s → 17 s), `evaluate_predictions()` sur 3,6e7 cellules de **2,18 Go à
  1,20 Go**.
* **`evaluate_predictions()` — `cor()` sur des vecteurs de longueurs
  différentes** — les valeurs de `prediction` et de `reference` étaient
  filtrées **séparément** par `values(x, na.rm = TRUE)` avant d'être passées à
  `cor()`. Dès que les deux masques de `NA` différaient — un CHM de référence
  troué, une prédiction rognée — les deux vecteurs n'avaient pas la même
  longueur et la fonction s'arrêtait sur une erreur. L'échantillonnage par
  paires écarte les cellules `NA` de part ou d'autre, donc le cas fonctionne.
* **`setup_conda_env()` ne vérifiait qu'un quart du stack Python** — le
  contrôle portait sur `torch, numpy, rasterio, huggingface_hub`, alors que
  `pipeline_aoi_to_chm()` charge aussi `geopandas, shapely, pyproj, rioxarray,
  xarray` — plus `torchvision, timm, segmentation_models_pytorch` pour
  l'inférence. Un environnement créé à la main avec seulement torch/rasterio
  passait le diagnostic, puis échouait sur `ModuleNotFoundError: geopandas` à
  l'étape AOI/CHM, sans produire de `chm_1_5m.tif`. La liste complète est
  désormais centralisée dans `OPEN_CANOPY_PY_MODULES` (nom d'import → nom pip)
  et vérifiée en entier.
* **CRS non rattachable sur tous les rasters produits** — les GeoTIFF renvoyés
  par le WMS IGN portent un WKT dont le *nom* est `"EPSG:2154"` mais qui n'a
  pas de bloc d'autorité `ID["EPSG",2154]` : le datum y est `"unnamed"` et
  l'ellipsoïde `"unretrievable - using WGS84"`. Le garde `is.na(crs(r)) ||
  crs(r) == ""` ne voyait pas passer ce cas, et le WKT dégénéré se propageait
  jusqu'à `chm_predicted_0_2m.tif`, où `sf::st_crs(x)$epsg` lit `NA`. Tout aval
  qui écrit un GeoPackage depuis ces rasters — segmentation de houppiers,
  exports vectoriels — embarquait un CRS non rattachable. `ancrer_crs_l93()`
  re-tamponne à `EPSG:2154` les seuls rasters dépourvus de code d'autorité (un
  raster correctement identifié dans un autre CRS est laissé tel quel) et est
  appliqué en quatre points : la tuile WMS — la racine du problème —, la
  branche cache des mosaïques (ce qui rattrape les `ortho_*.tif` déjà sur
  disque), leur écriture, et le CHM au retour de l'inférence Python, dont
  l'aller-retour par rasterio recopie le CRS d'entrée.

### Changements de comportement

* **Statistiques descriptives : exactes ou estimées, selon ce que `global()`
  sait faire.** Min, max, moyenne, écart-type, comptages, biais, RMSE et MAE
  restent **exacts**. En revanche les quantiles (médiane, Q25, Q75), les parts
  de surface (`pct_forest`, `pct_tall_trees`, `pct_vegetation`,
  `pct_dense_veg`, `pct_bare_soil`) et le R² sont désormais estimés sur un
  échantillon régulier borné à `max_cells` cellules — et redeviennent exacts
  dès que le raster tient dans cette limite. Écart mesuré sur les sorties IGN
  réelles à 0,20 m : **0,008 point de pourcentage** sur les parts, moins d'un
  millimètre sur les quantiles de hauteur, 1e-4 sur le R². La voie « exacte »
  évidente pour les parts, `global(chm >= 2, "mean")`, a été mesurée puis
  écartée : elle matérialise une couche de plus et coûte **6 Go, soit plus que
  `values()`** — c'est écrit dans la documentation des fonctions pour que le
  piège ne soit pas retenté.
* **Attention à `global(x, "rms")`** : terra y divise par *n−1*, pas par *n*.
  Le prendre pour la RMSE introduit un biais silencieux ; elle est reconstruite
  depuis les moments, `sqrt(sd_pop² + moyenne²)`.

### New Features

* **`setup_conda_env(install_missing = TRUE)`** installe via pip les modules
  Python manquants dans l'environnement. Par défaut (`FALSE`), la fonction se
  contente de les signaler avec la commande exacte à lancer, et retourne
  invisiblement le vecteur des modules absents.
* `ancrer_crs_l93()` est exportée : elle permet de réparer à la lecture les
  rasters produits par une version antérieure, sans relancer le pipeline.
* `compute_chm_stats()`, `compute_irc_stats()` et `evaluate_predictions()`
  acceptent `max_cells` (défaut `1e6`), et `plot_chm_histogram()` /
  `plot_prediction_comparison()` acceptent respectivement `max_cells` et
  `n_points`, pour régler le compromis précision/mémoire.

# opencanopy 0.1.3

Patch release ciblée sur l'empreinte mémoire des statistiques rasters.

### Bug Fixes

* **`pipeline_aoi_to_chm()` — OOM à la dernière étape sur les grandes AOI** —
  le pourcentage de végétation était calculé par
  `sum(values(clean_mask, na.rm = TRUE)) / sum(!is.na(values(clean_mask)))`,
  qui rapatriait le masque **entier** en mémoire, deux fois, plus la copie
  temporaire du `!is.na()`. Sur Couchey (535,6 ha, 28 481 × 14 695 = 418 528 295
  cellules à 0,20 m), cela représentait plusieurs Go de vecteurs R au pic de
  mémoire du pipeline — juste après un `mask()` qui, lui, streamait correctement
  vers le disque. Le processus était tué par l'OOM killer après 3 h 20 de CPU,
  tous les livrables écrits mais aucun indicateur calculé. Le calcul passe
  désormais par `global(clean_mask, "mean", na.rm = TRUE)`, qui streame par
  blocs : sur un masque logique, la moyenne en ignorant les `NA` **est** la
  proportion cherchée, au bit près.
* **Même motif ailleurs** — `values()` sur une couche pleine résolution pour un
  simple `min`/`max`/`mean`/`sum` a été remplacé par `global()` dans
  `pipeline_aoi_to_chm()` (statistiques du CHM, calculées une seule fois au lieu
  de quatre lectures complètes), `load_chm()`, `compute_ndvi()`,
  `compute_ndwi()`, `mask_vegetation()`, `compute_irc_stats()` (moyennes des
  trois bandes de l'ortho IRC), `detect_canopy_loss()` et `export_to_gpkg()`
  (test binaire). `compute_chm_stats()` déduit `n_na` de `ncell()` au lieu d'une
  seconde lecture complète. Résultats inchangés dans tous les cas, y compris sur
  les couches entièrement `NA`.

# opencanopy 0.1.2

Patch release ciblée sur la robustesse des téléchargements WMS IGN.

### Bug Fixes

* **WMS IGN — tuiles NULL après échec HTTP 502** — quand une tuile
  échouait (HTTP 5xx, timeout), `tile_rasters[[idx]] <- r` était
  sauté mais `idx` était quand même incrémenté, laissant des trous
  `NULL` dans la liste. `terra::merge` plantait alors avec
  `[sprc] list element is a: NULL`. `download_ign_tiled()` filtre
  désormais les `NULL` avant le mosaïquage et émet un avertissement
  clair indiquant combien de tuiles ont échoué.
* **WMS IGN — retry sur erreurs serveur transitoires** —
  `download_wms_tile()` retente jusqu'à 4 fois avec backoff
  exponentiel (1s, 2s, 4s, 8s) sur les erreurs HTTP 502/503/504,
  timeouts et resets. Élimine en pratique la plupart des échecs
  liés à la charge de la Géoplateforme IGN.

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
