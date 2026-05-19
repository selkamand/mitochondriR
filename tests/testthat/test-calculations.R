test_that("mtdna_per_cell calculation works", {
  expect_equal(
    calculate_mtdna(
      mitochondrial_depth = 100,
      autosome_depth = 10,
      ploidy = 2,
      per_cell = TRUE
    ),
    40
  )

  expect_equal(
    calculate_mtdna(
      mitochondrial_depth = 100,
      autosome_depth = 10,
      ploidy = 2,
      per_cell = FALSE
    ),
    20
  )


  expect_equal(
    calculate_mtdna(
      mitochondrial_depth = 100,
      autosome_depth = 10,
      ploidy = 2,
      per_cell = FALSE
    ),
    20
  )

})

test_that("mtdna_per_cell calculation works with tumour correction", {

  # Diploid tumour
  expect_equal(
    calculate_mtdna_tumour(
      mitochondrial_depth = 100,
      autosome_depth = 10,
      background_ploidy = 2,
      tumour_ploidy = 2,
      tumour_purity = 1,
      per_cell = TRUE
    ),
    40
  )

  # Genome-doubled tumour sample
  expect_equal(
    calculate_mtdna_tumour(
      mitochondrial_depth = 100,
      autosome_depth = 10,
      background_ploidy = 2,
      tumour_ploidy = 4,
      tumour_purity = 1,
      per_cell = TRUE
    ),
    160
  )

  # Ignore tumour ploidy if purity is 0
  expect_equal(
    calculate_mtdna_tumour(
      mitochondrial_depth = 100,
      autosome_depth = 10,
      background_ploidy = 2,
      tumour_ploidy = 100,
      tumour_purity = 0,
      per_cell = TRUE
    ),
    40
  )


  # 50% purity tetraploid tumour with 50% diplopid background.
  expect_equal(
    calculate_mtdna_tumour(
      mitochondrial_depth = 100,
      autosome_depth = 10,
      background_ploidy = 2,
      tumour_ploidy = 4,
      tumour_purity = 0.5,
      per_cell = TRUE
    ),
    90
  )
})


test_that("mtdna from idxstats works", {
  path <- system.file("idxstats", package = "mitochondriR")
  expect_equal(
    round(mtdna_from_idxstats(
      path = path,
      tumour_purity = 1,
      tumour_ploidy = 2,
      background_ploidy = 2,
      per_cell = TRUE,
      verbose = FALSE
    ), digits=0),
    47806
  )
})
