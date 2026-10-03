library(testthat)


test_that("make_mep creates deterministic epoched data with known geometry", {
  first <- make_mep(
    n_time = 600,
    n_channels = 2,
    n_trials = 4,
    sr = 5000,
    mep_latency_ms = 20,
    mep_dur_ms = 10,
    seed = 10
  )
  second <- make_mep(
    n_time = 600,
    n_channels = 2,
    n_trials = 4,
    sr = 5000,
    mep_latency_ms = 20,
    mep_dur_ms = 10,
    seed = 10
  )
  data <- SummarizedExperiment::assay(first)

  expect_equal(dim(data), c(600, 2, 4))
  expect_equal(data, SummarizedExperiment::assay(second))
  expect_equal(
    SummarizedExperiment::colData(first)$label,
    c("MEP1", "MEP2")
  )
  expect_error(
    make_mep(n_time = 200),
    "increase n_time"
  )
})


test_that("mepAmplitude recovers planted peak-to-peak amplitude", {
  experiment <- make_mep(
    n_trials = 3,
    sr = 5000,
    plateau = 2,
    s50 = 50,
    k = 5,
    intensities = rep(50, 3),
    baseline_sd = 0,
    seed = 3
  )
  amplitude <- mepAmplitude(
    experiment,
    response_window_ms = c(10, 60),
    epoch_start_ms = -20
  )

  expect_equal(amplitude$p2p, rep(1, 3), tolerance = 0.01)
  expect_equal(nrow(amplitude), 3L)
  expect_true(all(amplitude$pos_peak_sample %in% c(213L, 214L)))
  expect_true(all(amplitude$neg_peak_sample == 238L))
  expect_true(all(amplitude$pos_peak_time_ms %in% c(22.4, 22.6)))
  expect_equal(amplitude$intensity_pct_mso, rep(50, 3))
})


test_that("mepAmplitude promotes 2-D input and selects channels by label", {
  simulation <- make_mep(
    n_channels = 2,
    n_trials = 2,
    baseline_sd = 0,
    seed = 4
  )
  first_trial <- SummarizedExperiment::assay(simulation)[, , 1]
  experiment <- PhysioCore::PhysioExperiment(
    assays = list(raw = first_trial),
    colData = S4Vectors::DataFrame(
      label = c("FDI", "APB"),
      type = c("MEP", "MEP")
    ),
    samplingRate = 5000
  )
  experiment <- setTMSpulses(experiment, 40)

  amplitude <- mepAmplitude(
    experiment,
    epoch_start_ms = -20,
    channels = "APB"
  )
  expect_equal(nrow(amplitude), 1L)
  expect_equal(amplitude$channel, 2L)
  expect_identical(amplitude$channel_label, "APB")
  expect_equal(amplitude$trial, 1L)
})


test_that("mepAmplitude promotes a 2-D assay without dimnames", {
  simulation <- make_mep(
    n_channels = 1,
    n_trials = 1,
    baseline_sd = 0,
    seed = 5
  )
  values <- unname(SummarizedExperiment::assay(simulation)[, , 1])
  experiment <- PhysioCore::PhysioExperiment(
    assays = list(raw = matrix(values, ncol = 1)),
    colData = S4Vectors::DataFrame(label = "FDI", type = "MEP"),
    samplingRate = 5000
  )

  amplitude <- mepAmplitude(experiment, epoch_start_ms = -20)
  expect_equal(nrow(amplitude), 1L)
  expect_true(is.finite(amplitude$p2p))
})


test_that("MEP summaries return typed empty tables for zero dimensions", {
  zero_channel <- PhysioCore::PhysioExperiment(
    assays = list(raw = array(numeric(), c(10, 0, 1))),
    samplingRate = 1000
  )
  zero_trial <- PhysioCore::PhysioExperiment(
    assays = list(raw = array(numeric(), c(10, 1, 0))),
    colData = S4Vectors::DataFrame(label = "FDI", type = "MEP"),
    samplingRate = 1000
  )

  for (experiment in list(zero_channel, zero_trial)) {
    amplitude <- mepAmplitude(experiment)
    latency <- mepLatency(experiment)
    expect_equal(nrow(amplitude), 0L)
    expect_equal(nrow(latency), 0L)
    expect_type(amplitude$channel, "integer")
    expect_type(latency$latency_ms, "double")
  }
})


test_that("mepLatency is within one sample of planted onset", {
  sr <- 5000
  latency_ms <- 20
  epoch_start_ms <- -20
  experiment <- make_mep(
    n_trials = 8,
    sr = sr,
    mep_latency_ms = latency_ms,
    epoch_start_ms = epoch_start_ms,
    intensities = rep(65, 8),
    baseline_sd = 0.01,
    seed = 7
  )
  latency <- mepLatency(
    experiment,
    epoch_start_ms = epoch_start_ms,
    threshold_sd = 3
  )
  truth <- as.integer(round(
    (latency_ms - epoch_start_ms) / 1000 * sr
  )) + 1L

  expect_true(all(abs(latency$onset_sample - truth) <= 1L))
  expect_true(all(
    abs(latency$latency_ms - latency_ms) <= 1000 / sr + 1e-6
  ))
})


test_that("mepLatency handles flat and no-prestimulus epochs", {
  flat <- PhysioCore::PhysioExperiment(
    assays = list(raw = array(0, c(600, 1, 2))),
    colData = S4Vectors::DataFrame(label = "FDI", type = "MEP"),
    samplingRate = 5000
  )
  expect_warning(
    latency <- mepLatency(flat, epoch_start_ms = -20),
    "Degenerate baseline"
  )
  expect_true(all(is.na(latency$onset_sample)))
  expect_true(all(is.na(latency$latency_ms)))

  no_pre <- make_mep(
    n_time = 400,
    n_trials = 2,
    epoch_start_ms = 0,
    mep_latency_ms = 20,
    baseline_sd = 0.01,
    seed = 8
  )
  expect_warning(
    mepLatency(no_pre, epoch_start_ms = 0),
    "No pre-stimulus baseline"
  )
})


test_that("recruitmentCurve inverts a near-noiseless Boltzmann", {
  experiment <- make_mep(
    n_trials = 21,
    sr = 5000,
    plateau = 2,
    s50 = 50,
    k = 5,
    intensities = seq(30, 70, by = 2),
    baseline_sd = 1e-6,
    seed = 1
  )
  curve <- recruitmentCurve(
    experiment,
    channel = 1,
    epoch_start_ms = -20
  )

  expect_s3_class(curve, "recruitment_curve")
  expect_true(curve$converged)
  expect_equal(curve$plateau, 2, tolerance = 0.01)
  expect_equal(curve$s50, 50, tolerance = 0.01)
  expect_equal(curve$k, 5, tolerance = 0.01)
  expect_equal(curve$max_slope, 0.1, tolerance = 0.01)
  expect_gt(curve$r_squared, 0.999)
  expect_equal(curve$n, 21L)
  expect_equal(
    2 / (1 + exp((50 - 50) / 5)),
    1,
    tolerance = 1e-12
  )
})


test_that("recruitmentCurve recovers parameters within five percent", {
  experiment <- make_mep(
    n_trials = 25,
    sr = 5000,
    plateau = 2,
    s50 = 50,
    k = 5,
    intensities = seq(30, 70, length.out = 25),
    baseline_sd = 0.005,
    seed = 42
  )
  curve <- recruitmentCurve(experiment, epoch_start_ms = -20)

  expect_lt(abs(curve$plateau - 2) / 2, 0.05)
  expect_lt(abs(curve$s50 - 50) / 50, 0.05)
  expect_lt(abs(curve$k - 5) / 5, 0.05)
})


test_that("recruitmentCurve supports intensity overrides and guards inputs", {
  experiment <- make_mep(n_trials = 5, seed = 9)
  curve <- recruitmentCurve(
    experiment,
    epoch_start_ms = -20,
    intensity = c(10, 20, 30, 40, 50)
  )
  expect_equal(curve$data$intensity, c(10, 20, 30, 40, 50))
  expect_error(
    recruitmentCurve(
      experiment,
      epoch_start_ms = -20,
      intensity = rep(50, 5)
    ),
    "three distinct"
  )
  expect_error(
    recruitmentCurve(experiment, channel = 3, epoch_start_ms = -20),
    "invalid"
  )
})


test_that("MEP APIs validate windows and assay metadata", {
  experiment <- make_mep(n_trials = 3, seed = 11)

  expect_error(
    mepAmplitude(experiment, c(120, 130), epoch_start_ms = -20),
    "does not overlap"
  )
  expect_error(
    mepAmplitude(
      experiment,
      response_window_ms = c(60, 10),
      epoch_start_ms = -20
    ),
    "increasing"
  )
  expect_error(
    mepLatency(
      experiment,
      epoch_start_ms = -20,
      threshold_sd = 0
    ),
    "positive"
  )
  expect_error(
    mepAmplitude(
      experiment,
      epoch_start_ms = -20,
      channels = c(1, 1)
    ),
    "duplicate"
  )
})
