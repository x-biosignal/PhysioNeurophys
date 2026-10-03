# F-wave analysis

#' Empty per-trial F-wave table
#' @keywords internal
#' @noRd
.fwave_empty <- function() {
  data.frame(
    trial = integer(),
    detected = logical(),
    f_p2p = numeric(),
    f_latency_ms = numeric(),
    m_p2p = numeric()
  )
}


#' Summarize per-trial F-wave measurements
#' @keywords internal
#' @noRd
.fwave_summary <- function(waves) {
  n_trials <- nrow(waves)
  if (n_trials == 0L) {
    return(data.frame(
      n_trials = 0L,
      n_detected = 0L,
      persistence_pct = NA_real_,
      min_latency_ms = NA_real_,
      mean_latency_ms = NA_real_,
      max_latency_ms = NA_real_,
      chronodispersion_ms = NA_real_,
      mean_f_p2p = NA_real_,
      mmax = NA_real_,
      mean_f_m_ratio_pct = NA_real_
    ))
  }
  n_detected <- sum(waves$detected)
  detected_amplitude <- waves$f_p2p[waves$detected]
  detected_latency <- waves$f_latency_ms[
    waves$detected & is.finite(waves$f_latency_ms)
  ]
  mmax <- max(waves$m_p2p)
  if (length(detected_latency)) {
    min_latency <- min(detected_latency)
    mean_latency <- mean(detected_latency)
    max_latency <- max(detected_latency)
    chronodispersion <- max_latency - min_latency
  } else {
    min_latency <- mean_latency <- max_latency <- chronodispersion <- NA_real_
  }
  mean_f <- if (n_detected) mean(detected_amplitude) else NA_real_
  data.frame(
    n_trials = as.integer(n_trials),
    n_detected = as.integer(n_detected),
    persistence_pct = 100 * n_detected / n_trials,
    min_latency_ms = min_latency,
    mean_latency_ms = mean_latency,
    max_latency_ms = max_latency,
    chronodispersion_ms = chronodispersion,
    mean_f_p2p = mean_f,
    mmax = mmax,
    mean_f_m_ratio_pct = if (n_detected && mmax > 0) {
      100 * mean_f / mmax
    } else {
      NA_real_
    }
  )
}


#' Detect and summarize F-waves
#'
#' Measures late F-wave peak-to-peak amplitude and onset latency after
#' peripheral stimulation, then reports persistence, latency dispersion, and
#' the mean F/M amplitude ratio.
#'
#' @param x A `PhysioExperiment` with time x channels x trials data.
#' @param f_window_ms F-wave search window relative to the stimulus, in
#'   milliseconds.
#' @param m_window_ms Direct M-wave window relative to the stimulus, in
#'   milliseconds.
#' @param epoch_start_ms Time of the first epoch sample relative to the
#'   stimulus.
#' @param amplitude_threshold Minimum F-wave peak-to-peak amplitude in the
#'   assay's native unit.
#' @param baseline_window_ms Optional latency baseline interval. The default
#'   uses all pre-stimulus samples.
#' @param threshold_sd Baseline standard-deviation multiplier for onset.
#' @param min_consecutive_ms Required onset-threshold duration.
#' @param channel One channel label or 1-based index.
#' @param assay_name Optional assay name.
#'
#' @return An `f_wave_result` list containing per-trial `waves`, a one-row
#'   `summary`, and analysis `settings`.
#'
#' @references
#' Fisher MA (2007). F-waves: physiology and clinical uses.
#' *TheScientificWorldJournal*, 7:144-160. \doi{10.1100/tsw.2007.49}
#'
#' @export
#'
#' @examples
#' pe <- make_mep(n_trials = 5, mep_latency_ms = 25, seed = 3)
#' fWaveDetect(
#'   pe, m_window_ms = c(10, 18), f_window_ms = c(20, 60),
#'   epoch_start_ms = -20
#' )
fWaveDetect <- function(x,
                        f_window_ms = c(20, 60),
                        m_window_ms = c(2, 12),
                        epoch_start_ms = 0,
                        amplitude_threshold = 0.04,
                        baseline_window_ms = NULL,
                        threshold_sd = 3,
                        min_consecutive_ms = 0.5,
                        channel = 1,
                        assay_name = NULL) {
  source <- .mep_source(x, assay_name)
  epoch_start_ms <- .neuro_epoch_start(epoch_start_ms)
  if (!is.numeric(amplitude_threshold) ||
      length(amplitude_threshold) != 1L ||
      !is.finite(amplitude_threshold) || amplitude_threshold < 0) {
    stop("amplitude_threshold must be one non-negative finite value.",
         call. = FALSE)
  }
  if (!is.numeric(threshold_sd) || length(threshold_sd) != 1L ||
      !is.finite(threshold_sd) || threshold_sd <= 0 ||
      !is.numeric(min_consecutive_ms) ||
      length(min_consecutive_ms) != 1L ||
      !is.finite(min_consecutive_ms) || min_consecutive_ms <= 0) {
    stop("threshold_sd and min_consecutive_ms must be positive scalars.",
         call. = FALSE)
  }
  .neuro_baseline_window(baseline_window_ms)
  selected <- .neuro_channel(source, channel)
  windows <- .neuro_response_windows(
    m_window_ms,
    f_window_ms,
    epoch_start_ms,
    source,
    "m_window_ms",
    "f_window_ms"
  )
  settings <- list(
    f_window_ms = as.numeric(f_window_ms),
    m_window_ms = as.numeric(m_window_ms),
    epoch_start_ms = epoch_start_ms,
    amplitude_threshold = as.numeric(amplitude_threshold),
    baseline_window_ms = baseline_window_ms,
    threshold_sd = as.numeric(threshold_sd),
    min_consecutive_ms = as.numeric(min_consecutive_ms),
    channel = if (length(selected)) as.integer(selected) else integer(),
    assay_name = assay_name
  )
  if (!length(selected) || source$n_trials == 0L) {
    waves <- .fwave_empty()
    return(structure(
      list(
        waves = waves,
        summary = .fwave_summary(waves),
        settings = settings
      ),
      class = "f_wave_result"
    ))
  }

  latency <- mepLatency(
    x,
    response_window_ms = f_window_ms,
    epoch_start_ms = epoch_start_ms,
    baseline_window_ms = baseline_window_ms,
    threshold_sd = threshold_sd,
    min_consecutive_ms = min_consecutive_ms,
    channels = selected,
    assay_name = assay_name
  )
  rows <- vector("list", source$n_trials)
  for (trial in seq_len(source$n_trials)) {
    signal <- source$data[, selected, trial]
    m <- .mep_p2p(signal, windows$early)
    f <- .mep_p2p(signal, windows$late)
    detected <- unname(f["p2p"]) >= amplitude_threshold
    rows[[trial]] <- data.frame(
      trial = as.integer(trial),
      detected = detected,
      f_p2p = unname(f["p2p"]),
      f_latency_ms = if (detected) {
        latency$latency_ms[trial]
      } else {
        NA_real_
      },
      m_p2p = unname(m["p2p"])
    )
  }
  waves <- do.call(rbind, rows)
  structure(
    list(
      waves = waves,
      summary = .fwave_summary(waves),
      settings = settings
    ),
    class = "f_wave_result"
  )
}


#' @export
print.f_wave_result <- function(x, ...) {
  summary <- x$summary
  cat("F-wave analysis\n")
  cat(sprintf(
    "  Detected: %d/%d; persistence: %s%%\n",
    summary$n_detected,
    summary$n_trials,
    format(summary$persistence_pct, digits = 5)
  ))
  cat(sprintf(
    "  Mean latency: %s ms; chronodispersion: %s ms\n",
    format(summary$mean_latency_ms, digits = 5),
    format(summary$chronodispersion_ms, digits = 5)
  ))
  invisible(x)
}
