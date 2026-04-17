# Helper : construire un raster IRC synthétique (PIR, Rouge, Vert)
# conforme a la convention attendue par compute_ndvi/ndwi/gndvi/savi.

make_irc <- function(pir = 200, rouge = 80, vert = 100,
                     ncol = 10, nrow = 10) {
  b_pir   <- terra::rast(matrix(pir,   nrow = nrow, ncol = ncol))
  b_rouge <- terra::rast(matrix(rouge, nrow = nrow, ncol = ncol))
  b_vert  <- terra::rast(matrix(vert,  nrow = nrow, ncol = ncol))
  r <- c(b_pir, b_rouge, b_vert)
  names(r) <- c("PIR", "Rouge", "Vert")
  r
}
