library(testthat)


tms_eeg_fixture <- function(data, sr = 1000, labels = NULL,
                            pulse_absent = FALSE) {
  channels <- dim(data)[2L]
  if (is.null(labels)) {
    labels <- paste0("E", seq_len(channels))
  }
  experiment <- PhysioCore::PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = labels,
      type = rep("EEG", channels)
    ),
    samplingRate = sr
  )
  if (pulse_absent) {
    md <- S4Vectors::metadata(experiment)
    md$tms_eeg <- list(pulse_absent = TRUE)
    S4Vectors::metadata(experiment) <- md
  }
  experiment
}


tms_eeg_forward <- function(labels = paste0("E", 1:4)) {
  leadfield <- rbind(
    E1 = c(1.0, 0.0, 0.5),
    E2 = c(0.0, 1.0, -0.5),
    E3 = c(-0.7, 0.2, 1.0),
    E4 = c(-0.3, -1.2, -1.0)
  )
  leadfield <- leadfield[labels, , drop = FALSE]
  list(
    leadfield = leadfield,
    method = "analytic_fixture",
    coordinate_provenance = list(frame = "fixture")
  )
}


tms_eeg_flat <- function(data) {
  matrix(
    aperm(data, c(2L, 1L, 3L)),
    nrow = dim(data)[2L]
  )
}


tms_eeg_sound_reference <- function(data, leadfield, lambda,
                                    iterations, tolerance) {
  channels <- nrow(leadfield)
  H <- diag(channels) - matrix(1 / channels, channels, channels)
  Lr <- H %*% leadfield
  Yr <- H %*% tms_eeg_flat(data)
  sigma <- rep(1, channels)
  history <- matrix(sigma, nrow = 1L)
  metric <- numeric()
  converged <- FALSE
  for (iteration in seq_len(iterations)) {
    updated <- numeric(channels)
    for (sensor in seq_len(channels)) {
      minus <- setdiff(seq_len(channels), sensor)
      W <- diag(1 / sigma[minus], length(minus))
      Lt <- W %*% Lr[minus, , drop = FALSE]
      yt <- W %*% Yr[minus, , drop = FALSE]
      regularization <- lambda * sum(Lt^2) / length(minus)
      source <- t(Lt) %*% solve(
        Lt %*% t(Lt) + diag(regularization, length(minus)),
        yt
      )
      residual <- Yr[sensor, ] -
        as.numeric(Lr[sensor, , drop = FALSE] %*% source)
      updated[sensor] <- sqrt(mean(residual^2))
    }
    relative <- abs(updated - sigma) / sigma
    metric[iteration] <- max(relative)
    sigma <- updated
    history <- rbind(history, sigma)
    if (all(relative <= tolerance)) {
      converged <- TRUE
      break
    }
  }
  W <- diag(1 / sigma, channels)
  Lt <- W %*% Lr
  rank_h <- channels - 1L
  regularization <- lambda * sum(Lt^2) / rank_h
  K <- Lr %*% t(Lt) %*% solve(
    Lt %*% t(Lt) + diag(regularization, channels)
  ) %*% W
  list(
    history = history,
    metric = metric,
    sigma = sigma,
    K = K,
    clean = K %*% Yr,
    converged = converged,
    regularization = regularization
  )
}


test_that("pulse windows use rounded closed sample mapping at four rates", {
  for (sr in c(500, 1000, 1450, 5000)) {
    n_time <- as.integer(sr * 0.2) + 1L
    values <- seq_len(n_time)
    data <- array(values, c(n_time, 1, 1))
    experiment <- tms_eeg_fixture(data, sr)
    cleaned <- removeTMSpulse(
      experiment,
      pulse_window_ms = c(-1.1, 2.4),
      method = "linear",
      baseline_window_ms = c(-40, -20),
      epoch_start_ms = -50
    )
    step <- tail(S4Vectors::metadata(cleaned)$tms_eeg$steps, 1L)[[1L]]
    expected <- round((c(-1.1, 2.4) + 50) / 1000 * sr) + 1L

    expect_equal(
      unname(unlist(step$sample_ranges[1, ])),
      as.integer(expected)
    )
    expect_true(all(step$diagnostics$altered_sample_mask[
      seq.int(expected[1L], expected[2L]), , ,
      drop = FALSE
    ]))
    expect_identical(
      SummarizedExperiment::assay(experiment, "raw"),
      SummarizedExperiment::assay(cleaned, "raw")
    )
  }
})


test_that("linear and cubic bridges match independent arithmetic", {
  sr <- 1000
  samples <- seq_len(301)
  signal <- 0.002 * samples^2 - 0.3 * samples + 4
  data <- array(signal, c(301, 1, 1))
  experiment <- tms_eeg_fixture(data, sr)
  first <- 101L
  last <- 105L
  left <- first - 1L
  right <- last + 1L

  linear <- removeTMSpulse(
    experiment,
    pulse_window_ms = c(0, 4),
    method = "linear",
    baseline_window_ms = c(-80, -20),
    epoch_start_ms = -100
  )
  u <- (seq.int(first, last) - left) / (right - left)
  expected_linear <- (1 - u) * signal[left] + u * signal[right]
  expect_equal(
    SummarizedExperiment::assay(linear, "tms_pulse_clean")[
      seq.int(first, last), 1, 1
    ],
    expected_linear,
    tolerance = 1e-14
  )

  cubic <- removeTMSpulse(
    experiment,
    pulse_window_ms = c(0, 4),
    method = "cubic",
    support_ms = 3,
    baseline_window_ms = c(-80, -20),
    epoch_start_ms = -100
  )
  left_support <- seq.int(left - 3L, left)
  right_support <- seq.int(right, right + 3L)
  slope <- function(y, index) {
    centered <- index - mean(index)
    sum(centered * (y[index] - mean(y[index]))) / sum(centered^2)
  }
  ml <- slope(signal, left_support)
  mr <- slope(signal, right_support)
  width <- right - left
  h00 <- 2 * u^3 - 3 * u^2 + 1
  h10 <- u^3 - 2 * u^2 + u
  h01 <- -2 * u^3 + 3 * u^2
  h11 <- u^3 - u^2
  expected_cubic <- h00 * signal[left] + h10 * width * ml +
    h01 * signal[right] + h11 * width * mr
  expect_equal(
    SummarizedExperiment::assay(cubic, "tms_pulse_clean")[
      seq.int(first, last), 1, 1
    ],
    expected_cubic,
    tolerance = 1e-13
  )
})


test_that("recharge ranges merge and blank values are exact", {
  data <- array(0, c(401, 2, 2))
  data[, 1, 1] <- 3
  data[, 2, 1] <- -2
  data[, 1, 2] <- 7
  data[, 2, 2] <- 5
  experiment <- tms_eeg_fixture(data)
  recharges <- rbind(c(2.4, 5), c(20, 22))
  cleaned <- removeTMSpulse(
    experiment,
    pulse_window_ms = c(0, 2),
    recharge_windows_ms = recharges,
    method = "blank",
    blank_value = "baseline",
    baseline_window_ms = c(-80, -20),
    epoch_start_ms = -100
  )
  step <- tail(S4Vectors::metadata(cleaned)$tms_eeg$steps, 1L)[[1L]]
  output <- SummarizedExperiment::assay(cleaned, "tms_pulse_clean")

  expect_equal(
    step$sample_ranges,
    data.frame(
      first_sample = c(101L, 121L),
      last_sample = c(106L, 123L)
    )
  )
  expect_true(all(output[101:106, 1, 1] == 3))
  expect_true(all(output[101:106, 2, 1] == -2))
  expect_true(all(output[101:106, 1, 2] == 7))
  expect_true(all(output[101:106, 2, 2] == 5))

  zeroed <- removeTMSpulse(
    experiment,
    pulse_window_ms = c(0, 2),
    method = "blank",
    blank_value = "zero",
    baseline_window_ms = c(-80, -20),
    epoch_start_ms = -100
  )
  expect_true(all(
    SummarizedExperiment::assay(zeroed, "tms_pulse_clean")[101:103, , ] == 0
  ))
})


test_that("pulse removal preserves 2-D shape and guards invalid settings", {
  experiment <- tms_eeg_fixture(matrix(seq_len(401), 401, 1))
  output <- removeTMSpulse(
    experiment,
    pulse_window_ms = c(10, 12),
    method = "linear",
    baseline_window_ms = c(1, 5),
    epoch_start_ms = 0
  )
  expect_equal(
    dim(SummarizedExperiment::assay(output, "tms_pulse_clean")),
    c(401L, 1L)
  )
  expect_error(
    removeTMSpulse(experiment, method = "l"),
    "exactly one"
  )
  expect_error(
    removeTMSpulse(
      experiment, pulse_window_ms = c("10", "12"),
      method = "linear", baseline_window_ms = c(1, 5)
    ),
    "numeric"
  )
  expect_silent(
    removeTMSpulse(
      experiment, pulse_window_ms = c(10, 12),
      method = "linear",
      baseline_window_ms = c(-500, -400),
      epoch_start_ms = 0
    )
  )
  expect_error(
    removeTMSpulse(
      experiment, pulse_window_ms = c(500, 510),
      baseline_window_ms = c(1, 5)
    ),
    "does not overlap"
  )
  expect_error(
    removeTMSpulse(
      experiment, pulse_window_ms = c(0, 2),
      method = "cubic", support_ms = 0,
      baseline_window_ms = c(10, 20)
    ),
    "support_ms"
  )
  expect_error(
    removeTMSpulse(
      experiment, pulse_window_ms = c(0, 2),
      method = "blank", blank_value = "baseline",
      baseline_window_ms = c(0, 2)
    ),
    "wholly outside"
  )
  expect_error(
    removeTMSpulse(
      experiment, pulse_window_ms = c(10, 12),
      method = "linear", baseline_window_ms = c(1, 5),
      output_assay = "raw"
    ),
    "already exists"
  )
})


test_that("TMS-EEG audit metadata is appended without destructive coercion", {
  experiment <- tms_eeg_fixture(matrix(seq_len(201), 201, 1))
  md <- S4Vectors::metadata(experiment)
  md$tms_eeg <- "malformed"
  S4Vectors::metadata(experiment) <- md

  expect_error(
    removeTMSpulse(
      experiment,
      pulse_window_ms = c(10, 12),
      method = "linear",
      baseline_window_ms = c(1, 5)
    ),
    "must be a list"
  )
  expect_identical(
    S4Vectors::metadata(experiment)$tms_eeg,
    "malformed"
  )
})


test_that("SSP-SIR equals the declared source-informed reconstruction", {
  set.seed(20)
  labels <- paste0("E", 1:4)
  time <- seq(-100, 120, by = 1)
  trials <- 3L
  sources <- array(
    rnorm(3 * length(time) * trials, sd = 0.1),
    c(3, length(time), trials)
  )
  leadfield <- tms_eeg_forward(labels)$leadfield
  artifact_topography <- c(1, -1, 0.8, -0.8)
  artifact_time <- exp(-((time - 25) / 12)^2) * 8
  data <- array(0, c(length(time), 4, trials))
  for (trial in seq_len(trials)) {
    sensor <- leadfield %*% sources[, , trial] +
      artifact_topography %o% artifact_time
    data[, , trial] <- t(sensor)
  }
  experiment <- tms_eeg_fixture(
    data, labels = labels, pulse_absent = TRUE
  )
  result <- tmsSSPSIR(
    experiment,
    tms_eeg_forward(labels),
    artifact_window_ms = c(5, 50),
    n_artifact = 1,
    lambda = 0.05,
    epoch_start_ms = -100
  )
  step <- tail(S4Vectors::metadata(result)$tms_eeg$steps, 1L)[[1L]]
  H <- step$lead_field$reference_matrix
  R <- step$lead_field$cleaning_matrix
  expected <- R %*% H %*% tms_eeg_flat(data)
  observed <- tms_eeg_flat(
    SummarizedExperiment::assay(result, "tms_sspsir")
  )

  expect_equal(observed, expected, tolerance = 1e-12)
  expect_equal(dim(observed), c(4L, length(time) * trials))
  expect_equal(step$parameters$n_artifact, 1L)
  expect_equal(
    step$diagnostics$projector,
    diag(4) - step$diagnostics$artifact_topographies %*%
      t(step$diagnostics$artifact_topographies),
    tolerance = 1e-13
  )
  expect_gt(
    max(abs(observed - step$diagnostics$plain_ssp)),
    1e-6
  )
  expect_identical(
    SummarizedExperiment::assay(result, "raw"),
    SummarizedExperiment::assay(experiment, "raw")
  )
})


test_that("SSP-SIR reorders lead-field rows by exact labels", {
  set.seed(21)
  labels <- paste0("E", 1:4)
  data <- array(rnorm(201 * 4 * 2), c(201, 4, 2))
  experiment <- tms_eeg_fixture(
    data, labels = labels, pulse_absent = TRUE
  )
  ordered <- tmsSSPSIR(
    experiment, tms_eeg_forward(labels),
    artifact_window_ms = c(5, 30),
    n_artifact = 1,
    epoch_start_ms = -100
  )
  permutation <- c(3, 1, 4, 2)
  reordered_forward <- tms_eeg_forward(labels[permutation])
  reordered <- tmsSSPSIR(
    experiment, reordered_forward,
    artifact_window_ms = c(5, 30),
    n_artifact = 1,
    epoch_start_ms = -100
  )
  expect_equal(
    SummarizedExperiment::assay(ordered, "tms_sspsir"),
    SummarizedExperiment::assay(reordered, "tms_sspsir"),
    tolerance = 1e-12
  )

  invalid <- tms_eeg_forward(labels)
  rownames(invalid$leadfield)[1L] <- "missing"
  expect_error(
    tmsSSPSIR(
      experiment, invalid,
      artifact_window_ms = c(5, 30),
      n_artifact = 1,
      epoch_start_ms = -100
    ),
    "exactly cover"
  )

  positions_forward <- tms_eeg_forward(labels)
  positions_forward$electrode_positions <- data.frame(
    label = labels,
    x = seq_along(labels),
    y = 0,
    z = 0
  )
  rownames(positions_forward$leadfield) <- NULL
  via_positions <- tmsSSPSIR(
    experiment, positions_forward,
    artifact_window_ms = c(5, 30),
    n_artifact = 1,
    epoch_start_ms = -100
  )
  expect_equal(
    SummarizedExperiment::assay(ordered, "tms_sspsir"),
    SummarizedExperiment::assay(via_positions, "tms_sspsir"),
    tolerance = 1e-12
  )
})


test_that("SOUND matches an independent simultaneous-update implementation", {
  set.seed(22)
  labels <- paste0("E", 1:4)
  leadfield <- tms_eeg_forward(labels)$leadfield
  source <- array(rnorm(3 * 151 * 2), c(3, 151, 2))
  noise_sd <- c(0.03, 0.08, 0.15, 0.04)
  data <- array(0, c(151, 4, 2))
  for (trial in 1:2) {
    data[, , trial] <- t(
      leadfield %*% source[, , trial] +
        matrix(
          rnorm(4 * 151, sd = rep(noise_sd, each = 151)),
          4, 151, byrow = TRUE
        )
    )
  }
  experiment <- tms_eeg_fixture(
    data, labels = labels, pulse_absent = TRUE
  )
  lambda <- 0.1
  iterations <- 4L
  tolerance <- 0
  reference <- tms_eeg_sound_reference(
    data, leadfield, lambda, iterations, tolerance
  )
  result <- soundClean(
    experiment, tms_eeg_forward(labels),
    lambda = lambda,
    iterations = iterations,
    tolerance = tolerance,
    epoch_start_ms = -50
  )
  step <- tail(S4Vectors::metadata(result)$tms_eeg$steps, 1L)[[1L]]

  expect_equal(
    unname(step$diagnostics$noise_history),
    unname(reference$history),
    tolerance = 1e-11
  )
  expect_equal(
    step$diagnostics$convergence_metric,
    reference$metric,
    tolerance = 1e-11
  )
  expect_equal(
    step$lead_field$cleaning_matrix,
    reference$K,
    tolerance = 1e-11
  )
  expect_equal(
    tms_eeg_flat(SummarizedExperiment::assay(result, "tms_sound")),
    reference$clean,
    tolerance = 1e-10
  )
  expect_equal(
    step$lead_field$regularization,
    reference$regularization,
    tolerance = 1e-12
  )
  expect_false(step$diagnostics$converged)
})


test_that("SOUND is invariant to a common channel permutation", {
  set.seed(23)
  labels <- paste0("E", 1:4)
  data <- array(rnorm(121 * 4 * 2), c(121, 4, 2))
  experiment <- tms_eeg_fixture(
    data, labels = labels, pulse_absent = TRUE
  )
  original <- soundClean(
    experiment, tms_eeg_forward(labels),
    iterations = 3, tolerance = 0, epoch_start_ms = -40
  )

  permutation <- c(3, 1, 4, 2)
  permuted_data <- data[, permutation, , drop = FALSE]
  permuted_labels <- labels[permutation]
  permuted <- soundClean(
    tms_eeg_fixture(
      permuted_data, labels = permuted_labels, pulse_absent = TRUE
    ),
    tms_eeg_forward(permuted_labels),
    iterations = 3, tolerance = 0, epoch_start_ms = -40
  )
  restored <- SummarizedExperiment::assay(
    permuted, "tms_sound"
  )[, order(permutation), , drop = FALSE]
  expect_equal(
    restored,
    SummarizedExperiment::assay(original, "tms_sound"),
    tolerance = 1e-10
  )
})


test_that("source cleaners require pulse evidence and valid lead fields", {
  data <- array(0, c(201, 4, 1))
  experiment <- tms_eeg_fixture(data)
  forward <- tms_eeg_forward()
  expect_error(
    tmsSSPSIR(
      experiment, forward,
      artifact_window_ms = c(5, 30),
      n_artifact = 1,
      epoch_start_ms = -100
    ),
    "Remove the TMS pulse"
  )
  expect_error(
    soundClean(
      experiment, forward,
      epoch_start_ms = -100
    ),
    "Remove the TMS pulse"
  )

  declared <- tms_eeg_fixture(data, pulse_absent = TRUE)
  deficient <- forward
  deficient$leadfield[,] <- rep(c(1, -1, 1, -1), ncol(deficient$leadfield))
  expect_error(
    soundClean(declared, deficient, epoch_start_ms = -100),
    "rank at least two"
  )
  expect_error(
    soundClean(
      declared, forward,
      analysis_window_ms = c(500, 600),
      epoch_start_ms = -100
    ),
    "does not overlap"
  )
  expect_error(
    soundClean(
      declared, forward,
      iterations = .Machine$integer.max + 1,
      epoch_start_ms = -100
    ),
    "iterations"
  )
})


test_that("pulse provenance follows derived cleaning assays", {
  set.seed(24)
  labels <- paste0("E", 1:4)
  data <- array(rnorm(201 * 4 * 2), c(201, 4, 2))
  experiment <- tms_eeg_fixture(data, labels = labels)
  pulse_clean <- removeTMSpulse(
    experiment,
    pulse_window_ms = c(0, 2),
    method = "blank",
    blank_value = "zero",
    baseline_window_ms = c(-80, -20),
    epoch_start_ms = -100
  )
  ssp_clean <- tmsSSPSIR(
    pulse_clean,
    tms_eeg_forward(labels),
    artifact_window_ms = c(5, 30),
    n_artifact = 1,
    epoch_start_ms = -100,
    assay_name = "tms_pulse_clean"
  )
  expect_silent(
    sound <- soundClean(
      ssp_clean,
      tms_eeg_forward(labels),
      iterations = 1,
      tolerance = 0,
      epoch_start_ms = -100,
      assay_name = "tms_sspsir"
    )
  )
  provenance <- tail(
    S4Vectors::metadata(sound)$tms_eeg$steps, 1L
  )[[1L]]$pulse_provenance
  expect_identical(provenance$step$method, "removeTMSpulse")
  expect_false(provenance$hardware_declared)
})


test_that("source cleaners restore 2-D assay dimensions", {
  set.seed(25)
  labels <- paste0("E", 1:4)
  data <- matrix(rnorm(201 * 4), 201, 4)
  experiment <- tms_eeg_fixture(
    data, labels = labels, pulse_absent = TRUE
  )
  ssp <- tmsSSPSIR(
    experiment,
    tms_eeg_forward(labels),
    artifact_window_ms = c(5, 30),
    n_artifact = 1,
    epoch_start_ms = -100
  )
  sound <- soundClean(
    experiment,
    tms_eeg_forward(labels),
    iterations = 1,
    tolerance = 0,
    epoch_start_ms = -100
  )

  expect_equal(
    dim(SummarizedExperiment::assay(ssp, "tms_sspsir")),
    c(201L, 4L)
  )
  expect_equal(
    dim(SummarizedExperiment::assay(sound, "tms_sound")),
    c(201L, 4L)
  )
})


test_that("TEP averages recover canonical signed peaks per channel", {
  sr <- 1000
  epoch_start <- -200
  times <- seq(epoch_start, 250, by = 1)
  trials <- 5L
  data <- array(0, c(length(times), 2, trials))
  offsets <- cbind(seq_len(trials), -2 * seq_len(trials))
  for (trial in seq_len(trials)) {
    data[, 1, trial] <- offsets[trial, 1]
    data[, 2, trial] <- offsets[trial, 2]
  }
  peaks <- c(N15 = 15, P30 = 30, N45 = 45,
             P60 = 60, N100 = 100, P180 = 180)
  amplitudes <- c(N15 = -1, P30 = 2, N45 = -3,
                  P60 = 4, N100 = -5, P180 = 6)
  for (component in names(peaks)) {
    sample <- match(peaks[component], times)
    for (trial in seq_len(trials)) {
      data[sample, 1, trial] <- offsets[trial, 1] +
        amplitudes[component]
      data[sample, 2, trial] <- offsets[trial, 2] +
        2 * amplitudes[component]
    }
  }
  experiment <- tms_eeg_fixture(
    data, labels = c("C3", "C4"), pulse_absent = TRUE
  )
  result <- tepAverage(
    experiment,
    epoch_start_ms = epoch_start,
    channels = c("C4", "C3")
  )

  expect_s3_class(result, "tms_tep")
  expect_identical(
    unique(result$components$component),
    names(peaks)
  )
  expect_identical(
    result$components$channel_label[1:2],
    c("C4", "C3")
  )
  channel_one <- result$components[
    result$components$channel_label == "C3", ]
  expect_equal(channel_one$latency_ms, unname(peaks))
  expect_equal(channel_one$amplitude, unname(amplitudes))
  expect_false(any(channel_one$boundary_peak))
  expect_true(all(result$waveform$n_trials == trials))
})


test_that("TEP baseline precedes averaging and trimmed means are sample-wise", {
  times <- seq(-100, 100, by = 1)
  data <- array(0, c(length(times), 1, 5))
  offsets <- c(-4, -2, 0, 2, 4)
  peak <- match(30, times)
  for (trial in 1:5) {
    data[, 1, trial] <- offsets[trial]
    data[peak, 1, trial] <- offsets[trial] +
      c(2, 2, 2, 2, 100)[trial]
  }
  experiment <- tms_eeg_fixture(
    data, labels = "Cz", pulse_absent = TRUE
  )
  windows <- matrix(c(20, 40), nrow = 1L,
                    dimnames = list("P30", c("start", "end")))
  mean_result <- tepAverage(
    experiment,
    components = "P30",
    windows_ms = windows,
    baseline_window_ms = c(-80, -20),
    epoch_start_ms = -100,
    average = "mean"
  )
  trimmed <- tepAverage(
    experiment,
    components = "P30",
    windows_ms = windows,
    baseline_window_ms = c(-80, -20),
    epoch_start_ms = -100,
    average = "trimmed",
    trim = 0.2
  )

  expect_equal(mean_result$components$amplitude, 21.6)
  expect_equal(trimmed$components$amplitude, 2)
  expect_equal(trimmed$components$latency_ms, 30)
  expect_equal(
    trimmed$waveform$amplitude[trimmed$waveform$time_ms == -50],
    0
  )
})


test_that("TEP ties use earliest latency and missing provenance warns once", {
  times <- seq(-100, 100, by = 1)
  data <- array(0, c(length(times), 1, 2))
  data[times %in% c(25, 30), 1, ] <- 3
  experiment <- tms_eeg_fixture(data, labels = "Cz")
  windows <- matrix(c(20, 40), nrow = 1L,
                    dimnames = list("P30", c("start", "end")))
  expect_warning(
    result <- tepAverage(
      experiment,
      components = "P30",
      windows_ms = windows,
      baseline_window_ms = c(-80, -20),
      epoch_start_ms = -100
    ),
    "provenance is missing"
  )
  expect_equal(result$components$latency_ms, 25)
  expect_error(
    tepAverage(
      experiment, components = "P",
      baseline_window_ms = c(-80, -20),
      epoch_start_ms = -100
    ),
    "canonical"
  )
  expect_error(
    tepAverage(
      experiment, components = "P30",
      windows_ms = rbind(P30 = c(20, 40), N15 = c(10, 20)),
      baseline_window_ms = c(-80, -20),
      epoch_start_ms = -100
    ),
    "exactly covering"
  )
})


test_that("reference assets and Python arithmetic fixture are valid", {
  asset_dir <- system.file(
    "extdata", package = "PhysioNeurophys",
    mustWork = TRUE
  )
  manifest <- readLines(
    file.path(asset_dir, "tms_eeg_reference.sha256"),
    warn = FALSE
  )
  valid <- vapply(manifest, function(line) {
    fields <- strsplit(line, "[[:space:]]+")[[1L]]
    observed <- strsplit(
      system2(
        "sha256sum",
        file.path(asset_dir, fields[2L]),
        stdout = TRUE
      ),
      "[[:space:]]+"
    )[[1L]][1L]
    identical(observed, fields[1L])
  }, logical(1))
  fixture <- readRDS(file.path(asset_dir, "tms_eeg_reference.rds"))
  cross <- read.csv(
    file.path(asset_dir, "tms_eeg_cross_language.csv"),
    stringsAsFactors = FALSE
  )

  expect_true(all(valid))
  expect_identical(fixture$schema_version, 1L)
  expect_false(fixture$external$TESA_waveforms_included)
  expect_equal(dim(fixture$sspsir$leadfield), c(4L, 3L))

  signal <- matrix(seq_len(301L)^2, ncol = 1L)
  result <- removeTMSpulse(
    tms_eeg_fixture(signal, sr = 1450, labels = "E1"),
    pulse_window_ms = c(-1.7, 2.4),
    method = "linear",
    baseline_window_ms = c(-80, -30),
    epoch_start_ms = -100
  )
  bridge <- cross[cross$fixture == "linear_bridge", ]
  expect_equal(
    SummarizedExperiment::assay(
      result, "tms_pulse_clean"
    )[as.integer(bridge$key), 1L],
    bridge$value,
    tolerance = 1e-12
  )
})
