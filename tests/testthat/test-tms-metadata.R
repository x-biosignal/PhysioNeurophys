library(testthat)


test_that("TMS pulse metadata round-trips and preserves other metadata", {
  experiment <- make_mep(
    n_trials = 4,
    intensities = c(30, 40, 50, 60),
    seed = 1
  )
  metadata <- S4Vectors::metadata(experiment)
  metadata$provenance_marker <- "keep"
  S4Vectors::metadata(experiment) <- metadata

  experiment <- setTMSpulses(
    experiment,
    intensity_pct_mso = c(30, 40, 50, 60),
    target = "M1_FDI",
    coil = c("A", "B", "A", "B"),
    rmt_pct_mso = 42,
    onset_sec = c(1, 2, 3, 4)
  )
  pulses <- getTMSpulses(experiment)

  expect_identical(
    names(pulses),
    c(
      "pulse", "intensity_pct_mso", "target", "coil",
      "rmt_pct_mso", "onset_sec"
    )
  )
  expect_equal(pulses$pulse, 1:4)
  expect_equal(pulses$intensity_pct_mso, c(30, 40, 50, 60))
  expect_true(all(pulses$target == "M1_FDI"))
  expect_equal(pulses$rmt_pct_mso, rep(42, 4))
  expect_identical(
    S4Vectors::metadata(experiment)$provenance_marker,
    "keep"
  )
})


test_that("getTMSpulses returns a typed empty schema", {
  experiment <- PhysioCore::PhysioExperiment(
    assays = list(raw = matrix(0, 20, 1)),
    colData = S4Vectors::DataFrame(label = "FDI", type = "MEP"),
    samplingRate = 1000
  )
  pulses <- getTMSpulses(experiment)

  expect_equal(nrow(pulses), 0L)
  expect_type(pulses$pulse, "integer")
  expect_type(pulses$intensity_pct_mso, "double")
  expect_type(pulses$target, "character")
})


test_that("setTMSpulses validates trial alignment and ranges", {
  experiment <- make_mep(n_trials = 3, seed = 2)

  expect_error(
    setTMSpulses(experiment, c(30, 40)),
    "one finite value"
  )
  expect_error(
    setTMSpulses(experiment, c(30, 40, 101)),
    "\\[0, 100\\]"
  )
  expect_error(
    setTMSpulses(experiment, c(30, 40, 50), target = c("a", "b")),
    "target"
  )
  expect_error(
    setTMSpulses(experiment, c(30, 40, 50), onset_sec = -1),
    "onset_sec"
  )
})


test_that("make_mep writes one TMS row per trial", {
  experiment <- make_mep(
    n_trials = 5,
    intensities = seq(30, 70, length.out = 5),
    seed = 3
  )
  pulses <- getTMSpulses(experiment)

  expect_equal(nrow(pulses), 5L)
  expect_equal(
    pulses$intensity_pct_mso,
    seq(30, 70, length.out = 5)
  )
})
