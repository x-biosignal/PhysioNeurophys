# Synthetic MEP data

#' Simulate epoched TMS motor-evoked potentials
#'
#' Generates biphasic MEP waveforms with known latency and a Boltzmann
#' amplitude-intensity relationship.
#'
#' @param n_time Number of samples per epoch.
#' @param n_channels Number of MEP channels.
#' @param n_trials Number of TMS trials.
#' @param sr Sampling rate in Hz.
#' @param epoch_start_ms Time of the first sample relative to the TMS pulse.
#' @param intensities Per-trial stimulation intensity in `%MSO`.
#' @param plateau,s50,k Boltzmann recruitment parameters.
#' @param mep_latency_ms Planted MEP onset latency.
#' @param mep_dur_ms Duration of the one-cycle biphasic waveform.
#' @param baseline_sd Gaussian baseline-noise standard deviation.
#' @param seed Optional random seed.
#'
#' @return An epoched `PhysioExperiment` with TMS pulse metadata.
#' @export
#'
#' @examples
#' pe <- make_mep(n_trials = 5, seed = 1)
#' dim(SummarizedExperiment::assay(pe))
#' getTMSpulses(pe)
make_mep <- function(n_time = 600,
                     n_channels = 1,
                     n_trials = 12,
                     sr = 5000,
                     epoch_start_ms = -20,
                     intensities = NULL,
                     plateau = 2.0,
                     s50 = 50,
                     k = 5,
                     mep_latency_ms = 20,
                     mep_dur_ms = 10,
                     baseline_sd = 0.02,
                     seed = NULL) {
  integer_scalar <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < 1L ||
        value != as.integer(value)) {
      stop(sprintf("%s must be a positive integer.", name),
           call. = FALSE)
    }
    as.integer(value)
  }
  finite_scalar <- function(value, name, positive = FALSE,
                            nonnegative = FALSE) {
    valid <- is.numeric(value) && length(value) == 1L &&
      is.finite(value)
    if (positive) {
      valid <- valid && value > 0
    }
    if (nonnegative) {
      valid <- valid && value >= 0
    }
    if (!valid) {
      stop(sprintf("%s has an invalid numeric value.", name),
           call. = FALSE)
    }
    as.numeric(value)
  }
  n_time <- integer_scalar(n_time, "n_time")
  n_channels <- integer_scalar(n_channels, "n_channels")
  n_trials <- integer_scalar(n_trials, "n_trials")
  sr <- finite_scalar(sr, "sr", positive = TRUE)
  epoch_start_ms <- finite_scalar(epoch_start_ms, "epoch_start_ms")
  plateau <- finite_scalar(plateau, "plateau", positive = TRUE)
  s50 <- finite_scalar(s50, "s50")
  k <- finite_scalar(k, "k", positive = TRUE)
  mep_latency_ms <- finite_scalar(
    mep_latency_ms, "mep_latency_ms", nonnegative = TRUE
  )
  mep_dur_ms <- finite_scalar(mep_dur_ms, "mep_dur_ms", positive = TRUE)
  baseline_sd <- finite_scalar(
    baseline_sd, "baseline_sd", nonnegative = TRUE
  )
  if (is.null(intensities)) {
    intensities <- seq(30, 70, length.out = n_trials)
  }
  if (!is.numeric(intensities) || length(intensities) != n_trials ||
      any(!is.finite(intensities)) ||
      any(intensities < 0 | intensities > 100)) {
    stop("intensities must contain one finite value in [0, 100] per trial.",
         call. = FALSE)
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L ||
        !is.finite(seed) || seed != as.integer(seed)) {
      stop("seed must be NULL or one finite integer.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  onset <- as.integer(round(
    (mep_latency_ms - epoch_start_ms) / 1000 * sr
  )) + 1L
  duration <- max(1L, as.integer(round(mep_dur_ms / 1000 * sr)))
  if (onset < 1L || onset + duration > n_time) {
    stop(
      "The planted MEP lies outside the epoch; increase n_time or adjust timing.",
      call. = FALSE
    )
  }
  index <- seq.int(onset, onset + duration)
  waveform <- sin(
    2 * pi * (index - onset) / duration
  )
  amplitudes <- plateau /
    (1 + exp((s50 - intensities) / k))
  data <- array(
    stats::rnorm(
      n_time * n_channels * n_trials,
      mean = 0,
      sd = baseline_sd
    ),
    dim = c(n_time, n_channels, n_trials)
  )
  for (channel in seq_len(n_channels)) {
    for (trial in seq_len(n_trials)) {
      data[index, channel, trial] <- data[index, channel, trial] +
        0.5 * amplitudes[trial] * waveform
    }
  }
  labels <- paste0("MEP", seq_len(n_channels))
  dimnames(data) <- list(
    sample = NULL,
    channel = labels,
    trial = as.character(seq_len(n_trials))
  )
  experiment <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = labels,
      type = rep("MEP", n_channels)
    ),
    samplingRate = sr
  )
  setTMSpulses(
    experiment,
    intensity_pct_mso = intensities
  )
}
