library(testthat)


paired_experiment <- function(data, sr = 1000, labels = "FDI") {
  n_channels <- dim(data)[2L]
  PhysioCore::PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = labels[seq_len(n_channels)],
      type = rep("EMG", n_channels)
    ),
    samplingRate = sr
  )
}


paired_plant_p2p <- function(signal, sample, amplitude) {
  signal[sample] <- amplitude / 2
  signal[sample + 1L] <- -amplitude / 2
  signal
}


test_that("pairedPulse reproduces exact per-ISI ratios", {
  test_amplitudes <- rep(1, 4)
  conditioned <- c(0.4, 0.4, 0.5, 0.5, 1.5, 1.5, 1.4, 1.4)
  amplitudes <- c(test_amplitudes, conditioned)
  data <- array(0, c(100, 1, length(amplitudes)))
  for (trial in seq_along(amplitudes)) {
    data[, 1, trial] <- paired_plant_p2p(
      data[, 1, trial], 36, amplitudes[trial]
    )
  }
  experiment <- paired_experiment(data)
  condition <- c(rep("test", 4), rep("conditioned", 8))
  isi <- c(rep(NA_real_, 4), 2, 2, 3, 3, 10, 10, 15, 15)

  paired <- pairedPulse(
    experiment,
    condition,
    isi,
    epoch_start_ms = -20
  )
  summary <- siciIcf(paired)

  expect_equal(paired$isi_ms, c(2, 3, 10, 15))
  expect_equal(paired$ratio, c(0.4, 0.5, 1.5, 1.4), tolerance = 1e-12)
  expect_equal(paired$inhibition_pct, c(60, 50, -50, -40))
  expect_equal(paired$facilitation_pct, c(-60, -50, 50, 40))
  expect_equal(paired$n_conditioned, rep(2L, 4))
  expect_equal(paired$n_test, rep(4L, 4))
  expect_equal(summary$sici_ratio, 0.45)
  expect_equal(summary$sici_inhibition_pct, 55)
  expect_equal(summary$n_sici, 4L)
  expect_equal(summary$icf_ratio, 1.45)
  expect_equal(summary$icf_facilitation_pct, 45)
  expect_equal(summary$n_icf, 4L)
})


test_that("siciIcf weights by conditioned trials and includes boundaries", {
  paired <- data.frame(
    isi_ms = c(1, 5, 7, 20, 30),
    n_conditioned = c(1L, 3L, 2L, 4L, 1L),
    ratio = c(0.2, 0.6, 1.2, 1.6, 2)
  )
  result <- siciIcf(paired)

  expect_equal(result$sici_ratio, (0.2 + 3 * 0.6) / 4)
  expect_equal(result$n_sici, 4L)
  expect_equal(result$icf_ratio, (2 * 1.2 + 4 * 1.6) / 6)
  expect_equal(result$n_icf, 6L)

  missing <- siciIcf(
    paired[paired$isi_ms <= 5, ],
    icf_range_ms = c(7, 20)
  )
  expect_true(is.na(missing$icf_ratio))
  expect_true(is.na(missing$icf_facilitation_pct))
  expect_equal(missing$n_icf, 0L)
})


test_that("pairedPulse supports median aggregation and channel labels", {
  amplitudes <- c(1, 1, 1, 10, 0.5, 0.5)
  data <- array(0, c(100, 2, length(amplitudes)))
  for (trial in seq_along(amplitudes)) {
    data[, 2, trial] <- paired_plant_p2p(
      data[, 2, trial], 36, amplitudes[trial]
    )
  }
  experiment <- paired_experiment(data, labels = c("APB", "FDI"))
  condition <- c(rep("test", 4), rep("conditioned", 2))
  isi <- c(rep(NA_real_, 4), 2, 2)

  by_label <- pairedPulse(
    experiment, condition, isi,
    epoch_start_ms = -20,
    channel = "FDI",
    aggregate = "median"
  )
  by_index <- pairedPulse(
    experiment, condition, isi,
    epoch_start_ms = -20,
    channel = 2,
    aggregate = "median"
  )

  expect_equal(by_label, by_index)
  expect_equal(by_label$test_amplitude, 1)
  expect_equal(by_label$ratio, 0.5)
})


test_that("corticalSilentPeriod recovers exact onset and recovery", {
  sr <- 1000
  epoch_start_ms <- -100
  n_time <- 501
  times <- epoch_start_ms + (seq_len(n_time) - 1) / sr * 1000
  data <- array(rep(c(-1, 1), length.out = n_time), c(n_time, 1, 1))
  data[times >= 80 & times < 230, 1, 1] <- 0
  experiment <- paired_experiment(data, sr)

  result <- corticalSilentPeriod(
    experiment,
    search_window_ms = c(50, 300),
    epoch_start_ms = epoch_start_ms
  )

  expect_equal(result$baseline_rectified, 1)
  expect_equal(result$onset_ms, 80)
  expect_equal(result$offset_ms, 230)
  expect_equal(result$duration_ms, 150)
  expect_equal(result$onset_sample, 181L)
  expect_equal(result$offset_sample, 331L)
})


test_that("corticalSilentPeriod is sample-exact at 2 kHz", {
  sr <- 2000
  epoch_start_ms <- -100
  n_time <- 1001
  times <- epoch_start_ms + (seq_len(n_time) - 1) / sr * 1000
  data <- array(rep(c(-2, 2), length.out = n_time), c(n_time, 1, 1))
  data[times >= 75 & times < 225, 1, 1] <- 0

  result <- corticalSilentPeriod(
    paired_experiment(data, sr),
    search_window_ms = c(50, 300),
    epoch_start_ms = epoch_start_ms
  )

  expect_equal(result$onset_ms, 75)
  expect_equal(result$offset_ms, 225)
  expect_equal(result$duration_ms, 150)
})


test_that("CSP ignores short quiet and short return runs", {
  sr <- 1000
  epoch_start_ms <- -100
  n_time <- 501
  times <- epoch_start_ms + (seq_len(n_time) - 1) / sr * 1000
  data <- array(rep(c(-1, 1), length.out = n_time), c(n_time, 1, 1))
  data[times == 60, 1, 1] <- 0
  data[times >= 80 & times < 230, 1, 1] <- 0
  data[times >= 150 & times < 153, 1, 1] <- 1

  result <- corticalSilentPeriod(
    paired_experiment(data, sr),
    search_window_ms = c(50, 300),
    epoch_start_ms = epoch_start_ms,
    min_silence_ms = 10,
    min_return_ms = 5
  )

  expect_equal(result$onset_ms, 80)
  expect_equal(result$offset_ms, 230)
  expect_equal(result$duration_ms, 150)
})


test_that("CSP preserves missing onset and incomplete recovery", {
  sr <- 1000
  epoch_start_ms <- -100
  n_time <- 501
  times <- epoch_start_ms + (seq_len(n_time) - 1) / sr * 1000
  tonic <- array(rep(c(-1, 1), length.out = n_time), c(n_time, 1, 1))

  no_silence <- corticalSilentPeriod(
    paired_experiment(tonic, sr),
    search_window_ms = c(50, 300),
    epoch_start_ms = epoch_start_ms
  )
  expect_true(is.na(no_silence$onset_sample))
  expect_true(is.na(no_silence$duration_ms))

  no_return <- tonic
  no_return[times >= 80, 1, 1] <- 0
  incomplete <- corticalSilentPeriod(
    paired_experiment(no_return, sr),
    search_window_ms = c(50, 300),
    epoch_start_ms = epoch_start_ms
  )
  expect_equal(incomplete$onset_ms, 80)
  expect_true(is.na(incomplete$offset_sample))
  expect_true(is.na(incomplete$duration_ms))

  boundary_return <- no_return
  boundary_return[times >= 298, 1, 1] <- 1
  boundary <- corticalSilentPeriod(
    paired_experiment(boundary_return, sr),
    search_window_ms = c(50, 300),
    epoch_start_ms = epoch_start_ms,
    min_return_ms = 5
  )
  expect_true(is.na(boundary$offset_sample))
  expect_true(is.na(boundary$duration_ms))
})


test_that("CSP accepts an explicit baseline when no prestimulus epoch exists", {
  sr <- 1000
  n_time <- 401
  times <- (seq_len(n_time) - 1) / sr * 1000
  data <- array(rep(c(-1, 1), length.out = n_time), c(n_time, 1, 1))
  data[times >= 80 & times < 180, 1, 1] <- 0
  experiment <- paired_experiment(data, sr)

  expect_error(
    corticalSilentPeriod(
      experiment,
      search_window_ms = c(50, 300)
    ),
    "at least two"
  )
  result <- corticalSilentPeriod(
    experiment,
    search_window_ms = c(50, 300),
    baseline_window_ms = c(0, 40)
  )
  expect_equal(result$duration_ms, 100)
})


test_that("CSP promotes a 2-D assay without dimnames", {
  sr <- 1000
  epoch_start_ms <- -100
  n_time <- 501
  times <- epoch_start_ms + (seq_len(n_time) - 1) / sr * 1000
  signal <- rep(c(-1, 1), length.out = n_time)
  signal[times >= 80 & times < 180] <- 0
  experiment <- PhysioCore::PhysioExperiment(
    assays = list(raw = unname(matrix(signal, ncol = 1))),
    colData = S4Vectors::DataFrame(label = "FDI", type = "EMG"),
    samplingRate = sr
  )

  result <- corticalSilentPeriod(
    experiment,
    search_window_ms = c(50, 300),
    epoch_start_ms = epoch_start_ms
  )
  expect_equal(nrow(result), 1L)
  expect_equal(result$duration_ms, 100)
})


test_that("paired-pulse and CSP inputs are guarded", {
  experiment <- paired_experiment(array(0, c(100, 1, 3)))

  expect_error(
    pairedPulse(
      experiment,
      c("test", "conditioned", "conditioned"),
      c(NA, 2, 2),
      epoch_start_ms = -20
    ),
    "test MEP amplitude"
  )
  expect_error(
    pairedPulse(experiment, rep("test", 3), rep(NA_real_, 3)),
    "Both test"
  )
  expect_error(
    pairedPulse(
      experiment,
      c("test", "conditioned", "conditioned"),
      c(NA, NA, 2)
    ),
    "positive and finite"
  )
  expect_error(
    pairedPulse(
      experiment,
      c("test", "conditioned", "conditioned"),
      c(NaN, 2, 3)
    ),
    "must be NA"
  )
  expect_error(
    pairedPulse(
      experiment,
      c("test", "conditioned", "bad"),
      c(NA, 2, 3)
    ),
    "condition"
  )
  expect_error(
    pairedPulse(
      experiment,
      c("test", "conditioned", "conditioned"),
      c(NA, 2, 3),
      response_window_ms = c(-10, 60)
    ),
    "post-stimulus"
  )

  overflow_data <- array(0, c(100, 1, 2))
  overflow_data[, 1, 1] <- paired_plant_p2p(
    overflow_data[, 1, 1], 36, 1e-320
  )
  overflow_data[, 1, 2] <- paired_plant_p2p(
    overflow_data[, 1, 2], 36, 1
  )
  expect_error(
    pairedPulse(
      paired_experiment(overflow_data),
      c("test", "conditioned"),
      c(NA, 2),
      epoch_start_ms = -20
    ),
    "ratio is not finite"
  )

  paired <- data.frame(
    isi_ms = c(2, 10),
    n_conditioned = c(2L, 2L),
    ratio = c(0.5, 1.5)
  )
  expect_error(siciIcf(transform(paired, ratio = -1)), "invalid")
  expect_error(
    siciIcf(transform(
      paired,
      n_conditioned = .Machine$integer.max + 1
    )),
    "invalid"
  )
  expect_error(
    siciIcf(paired, sici_range_ms = c(1, 10)),
    "must not overlap"
  )
  expect_error(siciIcf(rbind(paired, paired[1, ])), "one row")
})


test_that("CSP rejects flat baselines and malformed settings", {
  sr <- 1000
  epoch_start_ms <- -100
  data <- array(0, c(501, 1, 1))
  experiment <- paired_experiment(data, sr)

  expect_error(
    corticalSilentPeriod(experiment, epoch_start_ms = epoch_start_ms),
    "zero or non-finite"
  )
  expect_error(
    corticalSilentPeriod(
      experiment,
      epoch_start_ms = epoch_start_ms,
      threshold_fraction = 1
    ),
    "\\(0, 1\\)"
  )
  expect_error(
    corticalSilentPeriod(
      experiment,
      epoch_start_ms = epoch_start_ms,
      min_return_ms = 0
    ),
    "positive"
  )
  expect_error(
    corticalSilentPeriod(
      experiment,
      epoch_start_ms = epoch_start_ms,
      min_return_ms = .Machine$double.xmax
    ),
    "too large"
  )
  expect_error(
    corticalSilentPeriod(
      experiment,
      epoch_start_ms = epoch_start_ms,
      baseline_window_ms = c(60, 80)
    ),
    "must precede"
  )
})


test_that("CSP returns a typed empty table for zero trials", {
  experiment <- PhysioCore::PhysioExperiment(
    assays = list(raw = array(numeric(), c(501, 1, 0))),
    colData = S4Vectors::DataFrame(label = "FDI", type = "EMG"),
    samplingRate = 1000
  )
  result <- corticalSilentPeriod(experiment, epoch_start_ms = -100)

  expect_equal(nrow(result), 0L)
  expect_type(result$trial, "integer")
  expect_type(result$duration_ms, "double")
})
