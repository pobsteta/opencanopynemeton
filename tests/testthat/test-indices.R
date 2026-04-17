test_that("compute_ndvi renvoie la valeur attendue sur raster constant", {
  irc <- make_irc(pir = 200, rouge = 80)
  ndvi <- suppressMessages(compute_ndvi(irc))

  expect_s4_class(ndvi, "SpatRaster")
  expect_equal(names(ndvi), "NDVI")
  # NDVI attendu = (200 - 80) / (200 + 80) = 0.4285714
  expect_equal(unique(terra::values(ndvi))[1], (200 - 80) / (200 + 80),
               tolerance = 1e-6)
})

test_that("compute_ndvi est borne dans [-1, 1]", {
  irc <- make_irc(pir = 255, rouge = 0)   # NDVI extremum haut
  ndvi <- suppressMessages(compute_ndvi(irc))
  vals <- terra::values(ndvi)
  expect_true(all(vals >= -1 & vals <= 1, na.rm = TRUE))
  expect_equal(unique(vals)[1], 1.0)
})

test_that("compute_gndvi utilise bien la bande Vert", {
  irc <- make_irc(pir = 200, vert = 50)
  gndvi <- compute_gndvi(irc)

  expect_equal(names(gndvi), "GNDVI")
  expect_equal(unique(terra::values(gndvi))[1], (200 - 50) / (200 + 50),
               tolerance = 1e-6)
})

test_that("compute_savi correspond a la formule L=0.5", {
  irc <- make_irc(pir = 200, rouge = 80)
  savi <- compute_savi(irc, L = 0.5)

  expected <- ((200 - 80) / (200 + 80 + 0.5)) * (1 + 0.5)
  expect_equal(names(savi), "SAVI")
  expect_equal(unique(terra::values(savi))[1], expected, tolerance = 1e-6)
})

test_that("compute_ndwi detecte l'eau (PIR bas, Vert haut) et rejette la vegetation", {
  # Eau : NIR bas, Vert eleve -> NDWI positif
  irc_eau <- make_irc(pir = 30, vert = 120)
  ndwi_eau <- suppressMessages(compute_ndwi(irc_eau))
  expect_gt(unique(terra::values(ndwi_eau))[1], 0)

  # Vegetation : NIR eleve, Vert modere -> NDWI negatif
  irc_veg <- make_irc(pir = 200, vert = 80)
  ndwi_veg <- suppressMessages(compute_ndwi(irc_veg))
  expect_lt(unique(terra::values(ndwi_veg))[1], 0)

  expect_equal(names(ndwi_eau), "NDWI")
})

test_that("mask_vegetation produit un raster binaire", {
  irc <- make_irc(pir = 200, rouge = 80)   # NDVI ~ 0.43 > 0.3
  ndvi <- suppressMessages(compute_ndvi(irc))
  mask <- suppressMessages(mask_vegetation(ndvi, threshold = 0.3))

  expect_equal(names(mask), "vegetation")
  expect_true(all(unique(terra::values(mask)) %in% c(0, 1, TRUE, FALSE)))
  # Avec NDVI = 0.43 partout, tout doit etre vegetation
  expect_true(all(terra::values(mask) == 1))
})
