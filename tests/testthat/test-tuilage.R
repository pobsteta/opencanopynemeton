# Tuilage / recollement de l'inférence : garantir qu'aucune tuile ne manque
# et qu'il n'y a pas de couture aux limites (cf. effet de bord sur le CHM).

make_grille <- function(largeur_m, hauteur_m, res_m = 1.4, nlyrs = 4) {
  nc <- round(largeur_m / res_m)
  nr <- round(hauteur_m / res_m)
  r <- terra::rast(nrows = nr, ncols = nc,
                   xmin = 843650, xmax = 843650 + nc * res_m,
                   ymin = 6485750, ymax = 6485750 + nr * res_m,
                   nlyrs = nlyrs, crs = "EPSG:2154")
  terra::values(r) <- seq_len(terra::ncell(r) * nlyrs)
  r
}

decouper <- function(r, ...) {
  suppressMessages(opencanopy:::make_inference_tiles(r, ...))
}

test_that("toutes les tuiles ont les memes dimensions en pixels", {
  tuiles <- decouper(make_grille(4600, 2600))
  dims <- t(vapply(tuiles, function(x) c(terra::ncol(x), terra::nrow(x)),
                   numeric(2)))
  expect_equal(nrow(unique(dims)), 1)
})

test_that("le decoupage couvre l'emprise entiere, sans trou", {
  for (dims in list(c(400, 300), c(1000, 1000), c(1050, 900), c(2400, 1300))) {
    r <- make_grille(dims[1], dims[2])
    tuiles <- decouper(r)
    preds <- lapply(unname(tuiles), function(x) x[[1]])
    mos <- opencanopy:::mosaiquer_predictions(preds,
                                              attr(tuiles, "margin_x"),
                                              attr(tuiles, "margin_y"))
    expect_equal(as.vector(terra::ext(mos)), as.vector(terra::ext(r)))
    expect_equal(as.numeric(terra::global(is.na(mos), "sum")[1, 1]), 0)
  }
})

test_that("le fondu pondere reconstitue exactement une entree constante", {
  # Partition de l'unite : sum(w_i) = 1 partout, donc la moyenne ponderee
  # de tuiles issues du meme raster doit redonner ce raster.
  r <- make_grille(2400, 1300)
  tuiles <- decouper(r)
  preds <- lapply(unname(tuiles), function(x) x[[1]])
  mos <- opencanopy:::mosaiquer_predictions(preds,
                                            attr(tuiles, "margin_x"),
                                            attr(tuiles, "margin_y"))
  ecart <- terra::global(abs(mos - r[[1]]), "max", na.rm = TRUE)[1, 1]
  expect_lt(as.numeric(ecart), 1e-6)
})

test_that("une tuile entierement NA n'efface pas ses voisines", {
  r <- make_grille(1960, 1400)
  tuiles <- decouper(r)
  preds <- lapply(unname(tuiles), function(x) x[[1]])
  preds[[1]][] <- NA

  mos <- opencanopy:::mosaiquer_predictions(preds,
                                            attr(tuiles, "margin_x"),
                                            attr(tuiles, "margin_y"))
  n_na <- as.numeric(terra::global(is.na(mos), "sum")[1, 1])
  # Le recouvrement des tuiles voisines comble une partie du trou :
  # il reste strictement moins que la surface d'une tuile.
  expect_gt(n_na, 0)
  expect_lt(n_na, terra::ncell(preds[[1]]))
})

test_that(".tile_starts couvre toujours la derniere colonne", {
  for (n in c(10, 100, 714, 715, 1000, 3286)) {
    for (size in c(50, 714)) {
      if (n <= size) next
      s <- opencanopy:::.tile_starts(n, size, size - 20)
      expect_equal(s[1], 1L)
      expect_equal(s[length(s)], as.integer(n - size + 1L))
      expect_true(all(diff(s) > 0))
    }
  }
})
