library(testthat)


neurophys_fixture <- function(data, sr = 1000, labels = "Soleus") {
  if (length(dim(data)) == 2L) {
    n_channels <- ncol(data)
  } else {
    n_channels <- dim(data)[2L]
  }
  PhysioCore::PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = labels[seq_len(n_channels)],
      type = rep("EMG", n_channels)
    ),
    samplingRate = sr
  )
}


plant_biphasic <- function(signal, positive_sample, amplitude) {
  signal[positive_sample] <- amplitude / 2
  signal[positive_sample + 1L] <- -amplitude / 2
  signal
}


test_that("H-reflex recruitment recovers exact H and M amplitudes", {
  intensities <- c(10, 20, 30, 40, 50)
  h_truth <- c(0.1, 0.5, 1, 0.8, 0.2)
  m_truth <- c(0.2, 0.8, 1.5, 2, 2.5)
  data <- array(0, c(80, 1, length(intensities)))
  for (trial in seq_along(intensities)) {
    data[, 1, trial] <- plant_biphasic(data[, 1, trial], 6, m_truth[trial])
    data[, 1, trial] <- plant_biphasic(data[, 1, trial], 26, h_truth[trial])
  }
  experiment <- neurophys_fixture(data)

  recruitment <- hReflexRecruitment(
    experiment,
    intensities,
    m_window_ms = c(2, 12),
    h_window_ms = c(20, 45)
  )
  ratio <- hMaxMMax(recruitment)
  threshold <- hReflexThreshold(recruitment, criterion = 0.2)

  expect_equal(recruitment$m_p2p, m_truth, tolerance = 1e-12)
  expect_equal(recruitment$h_p2p, h_truth, tolerance = 1e-12)
  expect_equal(ratio$hmax, 1)
  expect_equal(ratio$mmax, 2.5)
  expect_equal(ratio$hmax_mmax_ratio, 0.4)
  expect_equal(ratio$hmax_intensity, 30)
  expect_equal(ratio$mmax_intensity, 50)
  expect_equal(threshold$cutoff, 0.5)
  expect_equal(threshold$threshold_intensity, 20)
})


test_that("H-reflex summaries handle repeats, ties, and zero Mmax", {
  recruitment <- data.frame(
    intensity = c(20, 10, 20, 30),
    h_p2p = c(0.4, 0.2, 0.8, 0.8),
    m_p2p = c(1, 0.5, 1, 1)
  )

  ratio <- hMaxMMax(recruitment)
  expect_equal(ratio$hmax_intensity, 20)
  expect_equal(ratio$mmax_intensity, 20)
  expect_equal(
    hReflexThreshold(recruitment, 0.6, "absolute")$threshold_intensity,
    20
  )

  zero_m <- transform(recruitment, m_p2p = 0)
  expect_true(is.na(hMaxMMax(zero_m)$hmax_mmax_ratio))
  expect_true(is.na(
    hReflexThreshold(zero_m, 0.05, "mmax")$threshold_intensity
  ))
})


test_that("H-reflex latency recovers stimulus-relative onset", {
  sr <- 1000
  epoch_start_ms <- -20
  data <- array(0, c(100, 1, 2))
  data[seq_len(20), 1, ] <- rep(c(-1e-4, 1e-4), 10)
  onset <- as.integer(round((25 - epoch_start_ms) / 1000 * sr)) + 1L
  data[onset:(onset + 2L), 1, ] <- 1
  experiment <- neurophys_fixture(data, sr)

  latency <- hReflexLatency(
    experiment,
    h_window_ms = c(20, 45),
    epoch_start_ms = epoch_start_ms,
    min_consecutive_ms = 1
  )

  expect_equal(latency$onset_sample, rep(onset, 2))
  expect_equal(latency$h_latency_ms, rep(25, 2))
})


test_that("F-wave fixture reproduces persistence and latency statistics", {
  sr <- 1000
  epoch_start_ms <- -20
  latencies <- c(25, 26, 24, 28, 27, 25, 29)
  data <- array(0, c(100, 1, 10))
  data[seq_len(20), 1, ] <- rep(c(-1e-4, 1e-4), 10)
  for (trial in seq_len(10)) {
    data[, 1, trial] <- plant_biphasic(data[, 1, trial], 26, 2)
  }
  for (trial in seq_along(latencies)) {
    onset <- as.integer(round(
      (latencies[trial] - epoch_start_ms) / 1000 * sr
    )) + 1L
    data[onset, 1, trial] <- 0.05
    data[onset + 1L, 1, trial] <- 0.05
    data[onset + 2L, 1, trial] <- -0.05
  }
  experiment <- neurophys_fixture(data, sr)

  result <- fWaveDetect(
    experiment,
    f_window_ms = c(20, 60),
    m_window_ms = c(2, 12),
    epoch_start_ms = epoch_start_ms,
    amplitude_threshold = 0.04,
    min_consecutive_ms = 1
  )
  summary <- result$summary

  expect_s3_class(result, "f_wave_result")
  expect_equal(summary$n_detected, 7L)
  expect_equal(summary$persistence_pct, 70)
  expect_equal(summary$min_latency_ms, 24)
  expect_equal(summary$mean_latency_ms, mean(latencies))
  expect_equal(summary$max_latency_ms, 29)
  expect_equal(summary$chronodispersion_ms, 5)
  expect_equal(summary$mean_f_p2p, 0.1)
  expect_equal(summary$mmax, 2)
  expect_equal(summary$mean_f_m_ratio_pct, 5)
  expect_true(all(is.na(result$waves$f_latency_ms[8:10])))
})


test_that("F-wave detection separates amplitude and onset criteria", {
  data <- array(0, c(100, 1, 2))
  data[seq_len(20), 1, ] <- rep(c(-1e-3, 1e-3), 10)
  data[26:27, 1, ] <- c(0.02, -0.02)
  data[46:48, 1, 1] <- c(0.03, 0.03, -0.03)
  experiment <- neurophys_fixture(data)

  result <- fWaveDetect(
    experiment,
    f_window_ms = c(20, 60),
    m_window_ms = c(2, 12),
    epoch_start_ms = -20,
    amplitude_threshold = 0.05,
    min_consecutive_ms = 1
  )

  expect_identical(result$waves$detected, c(TRUE, FALSE))
  expect_true(is.finite(result$waves$f_latency_ms[1]))
  expect_true(is.na(result$waves$f_latency_ms[2]))
  expect_equal(result$summary$persistence_pct, 50)
})


test_that("H-reflex and F-wave APIs preserve typed empty results", {
  zero_channel <- PhysioCore::PhysioExperiment(
    assays = list(raw = array(numeric(), c(50, 0, 1))),
    samplingRate = 1000
  )
  zero_trial <- PhysioCore::PhysioExperiment(
    assays = list(raw = array(numeric(), c(50, 1, 0))),
    colData = S4Vectors::DataFrame(label = "Soleus", type = "EMG"),
    samplingRate = 1000
  )

  recruitment <- hReflexRecruitment(zero_channel, 10)
  latency <- hReflexLatency(zero_channel)
  result <- fWaveDetect(zero_trial)

  expect_equal(nrow(recruitment), 0L)
  expect_type(recruitment$trial, "integer")
  expect_equal(nrow(latency), 0L)
  expect_type(latency$h_latency_ms, "double")
  expect_equal(nrow(result$waves), 0L)
  expect_type(result$waves$detected, "logical")
  expect_equal(result$summary$n_trials, 0L)
  expect_true(is.na(result$summary$persistence_pct))
  expect_error(
    hReflexRecruitment(
      zero_trial, numeric(),
      m_window_ms = c(20, 10)
    ),
    "increasing"
  )
  expect_error(
    hReflexLatency(zero_trial, threshold_sd = 0),
    "positive"
  )
  expect_error(
    fWaveDetect(zero_trial, f_window_ms = c(60, 20)),
    "increasing"
  )
  expect_error(
    hReflexLatency(zero_trial, baseline_window_ms = 1),
    "two increasing"
  )
  expect_error(
    fWaveDetect(zero_trial, baseline_window_ms = c(1, NA)),
    "two increasing"
  )
})


test_that("channel labels and 2-D assays use the established MEP contract", {
  data <- matrix(0, nrow = 80, ncol = 2)
  data[, 2] <- plant_biphasic(data[, 2], 6, 1)
  data[, 2] <- plant_biphasic(data[, 2], 26, 0.5)
  experiment <- neurophys_fixture(data, labels = c("TA", "Soleus"))

  by_label <- hReflexRecruitment(
    experiment, 20, channel = "Soleus"
  )
  by_index <- hReflexRecruitment(
    experiment, 20, channel = 2
  )
  expect_equal(by_label, by_index)
  expect_equal(by_label$m_p2p, 1)
  expect_equal(by_label$h_p2p, 0.5)
})


test_that("H-reflex and F-wave inputs are guarded", {
  experiment <- neurophys_fixture(array(0, c(80, 1, 2)))
  table <- data.frame(intensity = 1, h_p2p = 1, m_p2p = 1)

  expect_error(
    hReflexRecruitment(experiment, 1),
    "one finite"
  )
  expect_error(
    hReflexRecruitment(
      experiment, c(1, 2),
      m_window_ms = c(2, 20),
      h_window_ms = c(20, 45)
    ),
    "not overlap"
  )
  expect_error(
    hReflexRecruitment(
      experiment, c(1, 2),
      m_window_ms = c(100, 110)
    ),
    "does not overlap"
  )
  expect_error(
    hReflexRecruitment(
      experiment, c(1, 2),
      m_window_ms = c(-2, 12)
    ),
    "post-stimulus"
  )
  expect_error(
    hMaxMMax(transform(table, h_p2p = -1)),
    "non-negative"
  )
  expect_error(
    hReflexThreshold(table, 1.1, "mmax"),
    "\\[0, 1\\]"
  )
  expect_error(
    fWaveDetect(experiment, amplitude_threshold = -1),
    "non-negative"
  )
  expect_error(
    fWaveDetect(experiment, threshold_sd = 0),
    "positive"
  )
  expect_error(
    fWaveDetect(experiment, channel = c(1, 1)),
    "duplicate"
  )
})


test_that("F-wave print is compact and returns invisibly", {
  data <- array(0, c(80, 1, 1))
  data[seq_len(20), 1, 1] <- rep(c(-1e-4, 1e-4), 10)
  result <- fWaveDetect(
    neurophys_fixture(data),
    epoch_start_ms = -20,
    amplitude_threshold = 1
  )
  output <- capture.output(returned <- print(result))

  expect_match(output[1], "F-wave analysis")
  expect_lt(length(output), 5L)
  expect_identical(returned, result)
})
