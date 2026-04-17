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
