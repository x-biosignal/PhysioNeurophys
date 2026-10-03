# Paired-pulse TMS and cortical silent-period analysis

#' Validate a post-stimulus analysis window
#' @keywords internal
#' @noRd
.paired_poststim_window <- function(window_ms,
                                    epoch_start_ms,
                                    source,
                                    name) {
  if (is.numeric(window_ms) && length(window_ms) == 2L &&
      all(is.finite(window_ms))) {
    mapped <- (window_ms - epoch_start_ms) / 1000 * source$sr
    if (any(!is.finite(mapped)) ||
        any(abs(mapped) > .Machine$integer.max - 1L)) {
      stop(sprintf("%s is outside the supported sample range.", name),
           call. = FALSE)
    }
  }
  samples <- .mep_window(
    window_ms,
    epoch_start_ms,
    source$sr,
    source$n_time,
    name
  )
  if (window_ms[1L] < 0) {
    stop(sprintf("%s must be a post-stimulus window.", name),
         call. = FALSE)
  }
  samples
}


#' Convert a minimum duration to a safe sample count
#' @keywords internal
#' @noRd
.paired_duration_samples <- function(duration_ms, sr, name) {
  samples <- ceiling(duration_ms / 1000 * sr)
  if (!is.finite(samples) || samples > .Machine$integer.max) {
    stop(sprintf("%s is too large for this sampling rate.", name),
         call. = FALSE)
  }
  max(1L, as.integer(samples))
}


#' Summarize paired-pulse MEP responses
#'
#' Compares conditioned motor-evoked-potential amplitudes with the aggregate
#' unconditioned test response at each interstimulus interval (ISI).
#'
#' @param x An epoched `PhysioExperiment` containing MEP trials.
#' @param condition Per-trial labels, exactly `"test"` or `"conditioned"`.
#' @param isi_ms Per-trial interstimulus interval in milliseconds. Use `NA` for
#'   test trials and a positive finite value for conditioned trials.
#' @param response_window_ms MEP response window relative to the test pulse.
#' @param epoch_start_ms Time of the first epoch sample relative to the test
#'   pulse.
#' @param channel One channel label or 1-based index.
#' @param assay_name Optional assay name.
#' @param aggregate Whether to summarize amplitudes by the mean or median.
#'
#' @return A data frame with one row per conditioned ISI and conditioned/test
#'   amplitudes, their ratio, inhibition, and facilitation.
#'
#' @references
#' Kujirai T, Caramia MD, Rothwell JC, et al. (1993). Corticocortical
#' inhibition in human motor cortex. *Journal of Physiology*, 471:501-519.
#' \doi{10.1113/jphysiol.1993.sp019912}
#'
#' @export
#'
#' @examples
#' pe <- make_mep(
#'   n_trials = 6,
#'   intensities = rep(50, 6),
#'   baseline_sd = 0,
#'   seed = 1
#' )
#' pairedPulse(
#'   pe,
#'   condition = c("test", "test", rep("conditioned", 4)),
#'   isi_ms = c(NA, NA, 2, 2, 10, 10),
#'   epoch_start_ms = -20
#' )
pairedPulse <- function(x,
                        condition,
                        isi_ms,
                        response_window_ms = c(10, 60),
                        epoch_start_ms = 0,
                        channel = 1,
                        assay_name = NULL,
                        aggregate = c("mean", "median")) {
  source <- .mep_source(x, assay_name)
  epoch_start_ms <- .neuro_epoch_start(epoch_start_ms)
  selected <- .neuro_channel(source, channel)
  .paired_poststim_window(
    response_window_ms,
    epoch_start_ms,
    source,
    "response_window_ms"
  )
  aggregate <- match.arg(aggregate)
  if (!length(selected)) {
    stop("channel must select exactly one channel.", call. = FALSE)
  }
  if (source$n_trials == 0L) {
    stop("pairedPulse requires at least one trial.", call. = FALSE)
  }
  if (!is.character(condition) ||
      length(condition) != source$n_trials ||
      anyNA(condition) ||
      any(!condition %in% c("test", "conditioned"))) {
    stop(
      "condition must contain one \"test\" or \"conditioned\" label per trial.",
      call. = FALSE
    )
  }
  if (!is.numeric(isi_ms) || length(isi_ms) != source$n_trials) {
    stop("isi_ms must contain one numeric value per trial.", call. = FALSE)
  }
  is_test <- condition == "test"
  is_conditioned <- condition == "conditioned"
  if (any(!is.na(isi_ms[is_test])) || any(is.nan(isi_ms[is_test]))) {
    stop("isi_ms must be NA for test trials.", call. = FALSE)
  }
  if (any(!is.finite(isi_ms[is_conditioned]) |
          isi_ms[is_conditioned] <= 0)) {
    stop("isi_ms must be positive and finite for conditioned trials.",
         call. = FALSE)
  }
  if (!any(is_test) || !any(is_conditioned)) {
    stop("Both test and conditioned trials are required.", call. = FALSE)
  }

  amplitudes <- mepAmplitude(
    x,
    response_window_ms = response_window_ms,
    epoch_start_ms = epoch_start_ms,
    channels = selected,
    assay_name = assay_name
  )$p2p
  if (any(!is.finite(amplitudes))) {
    stop("Measured MEP amplitudes must be finite.", call. = FALSE)
  }
  summarize <- if (aggregate == "mean") {
    mean
  } else {
    stats::median
  }
  test_amplitude <- summarize(amplitudes[is_test])
  if (!is.finite(test_amplitude) || test_amplitude <= 0) {
    stop("The aggregate test MEP amplitude must be positive and finite.",
         call. = FALSE)
  }

  levels <- sort(unique(isi_ms[is_conditioned]))
  rows <- vector("list", length(levels))
  for (index in seq_along(levels)) {
    at_level <- is_conditioned & isi_ms == levels[index]
    conditioned_amplitude <- summarize(amplitudes[at_level])
    ratio <- conditioned_amplitude / test_amplitude
    if (!is.finite(ratio)) {
      stop("The conditioned/test MEP ratio is not finite.", call. = FALSE)
    }
    rows[[index]] <- data.frame(
      isi_ms = as.numeric(levels[index]),
      n_conditioned = as.integer(sum(at_level)),
      conditioned_amplitude = conditioned_amplitude,
      n_test = as.integer(sum(is_test)),
      test_amplitude = test_amplitude,
      ratio = ratio,
      inhibition_pct = 100 * (1 - ratio),
      facilitation_pct = 100 * (ratio - 1)
    )
  }
  do.call(rbind, rows)
}


#' Summarize short-interval inhibition and intracortical facilitation
#'
#' @param paired A table returned by [pairedPulse()] containing `isi_ms`,
#'   `n_conditioned`, and `ratio`.
#' @param sici_range_ms Closed SICI ISI range in milliseconds.
#' @param icf_range_ms Closed ICF ISI range in milliseconds.
#'
#' @return A one-row data frame with conditioned-trial-weighted SICI and ICF
#'   ratios, percentages, and trial counts.
#' @export
#'
#' @examples
#' paired <- data.frame(
#'   isi_ms = c(2, 3, 10, 15),
#'   n_conditioned = rep(2L, 4),
#'   ratio = c(0.4, 0.5, 1.5, 1.4)
#' )
#' siciIcf(paired)
siciIcf <- function(paired,
                     sici_range_ms = c(1, 5),
                     icf_range_ms = c(7, 20)) {
  required <- c("isi_ms", "n_conditioned", "ratio")
  if (!is.data.frame(paired) || !all(required %in% names(paired))) {
    stop("paired must contain isi_ms, n_conditioned, and ratio columns.",
         call. = FALSE)
  }
  if (!is.numeric(paired$isi_ms) || any(!is.finite(paired$isi_ms)) ||
      any(paired$isi_ms <= 0) ||
      !is.numeric(paired$n_conditioned) ||
      any(!is.finite(paired$n_conditioned)) ||
      any(paired$n_conditioned <= 0) ||
      any(paired$n_conditioned != floor(paired$n_conditioned)) ||
      any(paired$n_conditioned > .Machine$integer.max) ||
      sum(paired$n_conditioned) > .Machine$integer.max ||
      !is.numeric(paired$ratio) || any(!is.finite(paired$ratio)) ||
      any(paired$ratio < 0)) {
    stop("paired contains invalid ISIs, counts, or ratios.", call. = FALSE)
  }
  if (anyDuplicated(paired$isi_ms)) {
    stop("paired must contain exactly one row per ISI.", call. = FALSE)
  }
  validate_range <- function(value, name) {
    if (!is.numeric(value) || length(value) != 2L ||
        any(!is.finite(value)) || value[1L] < 0 ||
        value[2L] <= value[1L]) {
      stop(sprintf("%s must be a non-negative increasing finite range.", name),
           call. = FALSE)
    }
    as.numeric(value)
  }
  sici_range_ms <- validate_range(sici_range_ms, "sici_range_ms")
  icf_range_ms <- validate_range(icf_range_ms, "icf_range_ms")
  if (!(sici_range_ms[2L] < icf_range_ms[1L] ||
        icf_range_ms[2L] < sici_range_ms[1L])) {
    stop("SICI and ICF ranges must not overlap.", call. = FALSE)
  }

  summarize <- function(range) {
    inside <- paired$isi_ms >= range[1L] &
      paired$isi_ms <= range[2L]
    count <- sum(paired$n_conditioned[inside])
    ratio <- if (count > 0) {
      sum(paired$ratio[inside] * paired$n_conditioned[inside]) / count
    } else {
      NA_real_
    }
    list(ratio = ratio, count = as.integer(count))
  }
  sici <- summarize(sici_range_ms)
  icf <- summarize(icf_range_ms)
  data.frame(
    sici_ratio = sici$ratio,
    sici_inhibition_pct = if (is.finite(sici$ratio)) {
      100 * (1 - sici$ratio)
    } else {
      NA_real_
    },
    n_sici = sici$count,
    icf_ratio = icf$ratio,
    icf_facilitation_pct = if (is.finite(icf$ratio)) {
      100 * (icf$ratio - 1)
    } else {
      NA_real_
    },
    n_icf = icf$count
  )
}


#' Empty cortical silent-period table
#' @keywords internal
#' @noRd
.csp_empty <- function() {
  data.frame(
    trial = integer(),
    channel = integer(),
    channel_label = character(),
    baseline_rectified = numeric(),
    onset_sample = integer(),
    offset_sample = integer(),
    onset_ms = numeric(),
    offset_ms = numeric(),
    duration_ms = numeric()
  )
}


#' Measure the cortical silent period
#'
#' Detects sustained suppression and recovery of rectified tonic EMG after a
#' TMS pulse using a baseline-relative threshold.
#'
#' @param x An epoched tonic-EMG `PhysioExperiment`.
#' @param search_window_ms Post-stimulus CSP search window.
#' @param epoch_start_ms Time of the first epoch sample relative to the pulse.
#' @param baseline_window_ms Optional quiet-reference interval. The default
#'   uses all pre-stimulus samples.
#' @param threshold_fraction Fraction of median rectified baseline below which
#'   EMG is considered silent.
#' @param min_silence_ms Required duration of suppression.
#' @param min_return_ms Required duration of recovered activity.
#' @param channel One channel label or 1-based index.
#' @param assay_name Optional assay name.
#'
#' @return A data frame with onset, offset, and duration for each trial.
#'
#' @references
#' Rossini PM, Burke D, Chen R, et al. (2015). Non-invasive electrical and
#' magnetic stimulation of the brain, spinal cord, roots and peripheral
#' nerves: basic principles and procedures for routine clinical and research
#' application. *Clinical Neurophysiology*, 126:1071-1107.
#' \doi{10.1016/j.clinph.2015.02.001}
#'
#' @export
#'
#' @examples
#' data <- array(rep(c(-1, 1), 250), c(500, 1, 1))
#' data[181:330, 1, 1] <- 0
#' pe <- PhysioExperiment(
#'   assays = list(raw = data),
#'   samplingRate = 1000
#' )
#' corticalSilentPeriod(pe, epoch_start_ms = -100)
corticalSilentPeriod <- function(x,
                                 search_window_ms = c(50, 400),
                                 epoch_start_ms = 0,
                                 baseline_window_ms = NULL,
                                 threshold_fraction = 0.25,
                                 min_silence_ms = 10,
                                 min_return_ms = 5,
                                 channel = 1,
                                 assay_name = NULL) {
  source <- .mep_source(x, assay_name)
  epoch_start_ms <- .neuro_epoch_start(epoch_start_ms)
  .neuro_baseline_window(baseline_window_ms)
  if (!is.numeric(threshold_fraction) ||
      length(threshold_fraction) != 1L ||
      !is.finite(threshold_fraction) ||
      threshold_fraction <= 0 || threshold_fraction >= 1) {
    stop("threshold_fraction must be one finite value in (0, 1).",
         call. = FALSE)
  }
  if (!is.numeric(min_silence_ms) || length(min_silence_ms) != 1L ||
      !is.finite(min_silence_ms) || min_silence_ms <= 0 ||
      !is.numeric(min_return_ms) || length(min_return_ms) != 1L ||
      !is.finite(min_return_ms) || min_return_ms <= 0) {
    stop("min_silence_ms and min_return_ms must be positive scalars.",
         call. = FALSE)
  }
  selected <- .neuro_channel(source, channel)
  search <- .paired_poststim_window(
    search_window_ms,
    epoch_start_ms,
    source,
    "search_window_ms"
  )
  times <- epoch_start_ms +
    (seq_len(source$n_time) - 1) / source$sr * 1000
  if (is.null(baseline_window_ms)) {
    baseline <- which(times < 0)
  } else {
    baseline <- which(
      times >= baseline_window_ms[1L] &
        times < baseline_window_ms[2L]
    )
  }
  if (length(baseline) < 2L) {
    stop("baseline_window_ms must select at least two samples.",
         call. = FALSE)
  }
  if (max(baseline) >= min(search)) {
    stop("The baseline interval must precede search_window_ms.",
         call. = FALSE)
  }
  if (!length(selected) || source$n_trials == 0L) {
    return(.csp_empty())
  }

  min_silence <- .paired_duration_samples(
    min_silence_ms,
    source$sr,
    "min_silence_ms"
  )
  min_return <- .paired_duration_samples(
    min_return_ms,
    source$sr,
    "min_return_ms"
  )
  rows <- vector("list", source$n_trials)
  for (trial in seq_len(source$n_trials)) {
    signal <- source$data[, selected, trial]
    baseline_rectified <- stats::median(abs(signal[baseline]))
    if (!is.finite(baseline_rectified) || baseline_rectified <= 0) {
      stop(
        sprintf("Trial %d has a zero or non-finite rectified baseline.", trial),
        call. = FALSE
      )
    }
    quiet <- abs(signal[search]) <=
      threshold_fraction * baseline_rectified
    onset <- .mep_onset(quiet, search, min_silence)
    if (is.na(onset)) {
      offset <- NA_integer_
    } else {
      after <- search[search > onset]
      offset <- .mep_onset(!quiet[search > onset], after, min_return)
    }
    onset_ms <- if (is.na(onset)) {
      NA_real_
    } else {
      epoch_start_ms + (onset - 1L) / source$sr * 1000
    }
    offset_ms <- if (is.na(offset)) {
      NA_real_
    } else {
      epoch_start_ms + (offset - 1L) / source$sr * 1000
    }
    rows[[trial]] <- data.frame(
      trial = as.integer(trial),
      channel = as.integer(selected),
      channel_label = source$labels[selected],
      baseline_rectified = baseline_rectified,
      onset_sample = onset,
      offset_sample = offset,
      onset_ms = onset_ms,
      offset_ms = offset_ms,
      duration_ms = if (is.na(onset) || is.na(offset)) {
        NA_real_
      } else {
        (offset - onset) / source$sr * 1000
      }
    )
  }
  do.call(rbind, rows)
}
