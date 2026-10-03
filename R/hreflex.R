# H-reflex analysis

#' Validate epoch timing
#' @keywords internal
#' @noRd
.neuro_epoch_start <- function(epoch_start_ms) {
  if (!is.numeric(epoch_start_ms) || length(epoch_start_ms) != 1L ||
      !is.finite(epoch_start_ms)) {
    stop("epoch_start_ms must be a finite numeric scalar.", call. = FALSE)
  }
  as.numeric(epoch_start_ms)
}


#' Resolve exactly one analysis channel
#' @keywords internal
#' @noRd
.neuro_channel <- function(source, channel) {
  if (source$n_channels == 0L) {
    return(integer())
  }
  selected <- .mep_channels(
    channel, source$labels, source$n_channels
  )
  if (length(selected) != 1L) {
    stop("channel must select exactly one channel.", call. = FALSE)
  }
  selected
}


#' Resolve ordered non-overlapping response windows
#' @keywords internal
#' @noRd
.neuro_response_windows <- function(early_window_ms,
                                    late_window_ms,
                                    epoch_start_ms,
                                    source,
                                    early_name,
                                    late_name) {
  early <- .mep_window(
    early_window_ms,
    epoch_start_ms,
    source$sr,
    source$n_time,
    early_name
  )
  late <- .mep_window(
    late_window_ms,
    epoch_start_ms,
    source$sr,
    source$n_time,
    late_name
  )
  if (early_window_ms[1L] < 0 || late_window_ms[1L] < 0) {
    stop(
      sprintf("%s and %s must be post-stimulus windows.",
              early_name, late_name),
      call. = FALSE
    )
  }
  if (max(early) >= min(late)) {
    stop(
      sprintf("%s must precede and not overlap %s.", early_name, late_name),
      call. = FALSE
    )
  }
  list(early = early, late = late)
}


#' Validate an optional baseline interval
#' @keywords internal
#' @noRd
.neuro_baseline_window <- function(baseline_window_ms) {
  if (!is.null(baseline_window_ms) &&
      (!is.numeric(baseline_window_ms) ||
       length(baseline_window_ms) != 2L ||
       any(!is.finite(baseline_window_ms)) ||
       baseline_window_ms[2L] <= baseline_window_ms[1L])) {
    stop("baseline_window_ms must contain two increasing finite times.",
         call. = FALSE)
  }
  baseline_window_ms
}


#' Empty H-reflex recruitment table
#' @keywords internal
#' @noRd
.hreflex_empty <- function() {
  data.frame(
    trial = integer(),
    intensity = numeric(),
    m_p2p = numeric(),
    m_pos_peak_sample = integer(),
    m_neg_peak_sample = integer(),
    h_p2p = numeric(),
    h_pos_peak_sample = integer(),
    h_neg_peak_sample = integer()
  )
}


#' Validate an H-reflex recruitment table
#' @keywords internal
#' @noRd
.hreflex_table <- function(recruitment) {
  required <- c("intensity", "h_p2p", "m_p2p")
  if (!is.data.frame(recruitment) ||
      !all(required %in% names(recruitment))) {
    stop(
      "recruitment must contain intensity, h_p2p, and m_p2p columns.",
      call. = FALSE
    )
  }
  for (name in required) {
    if (!is.numeric(recruitment[[name]]) ||
        any(!is.finite(recruitment[[name]]))) {
      stop(sprintf("%s must contain finite numeric values.", name),
           call. = FALSE)
    }
  }
  if (any(recruitment$h_p2p < 0 | recruitment$m_p2p < 0)) {
    stop("H- and M-wave amplitudes must be non-negative.", call. = FALSE)
  }
  recruitment
}


#' Measure an H-reflex recruitment series
#'
#' Measures peak-to-peak direct motor (M) and Hoffmann-reflex (H) responses in
#' separate post-stimulus windows for each trial.
#'
#' @param x A `PhysioExperiment` with time x channels x trials data.
#' @param intensity Finite stimulation intensity, one value per trial. Values
#'   remain in the caller's native stimulator unit.
#' @param m_window_ms Direct M-wave window relative to the stimulus, in
#'   milliseconds.
#' @param h_window_ms H-reflex window relative to the stimulus, in
#'   milliseconds.
#' @param epoch_start_ms Time of the first epoch sample relative to the
#'   stimulus.
#' @param channel One channel label or 1-based index.
#' @param assay_name Optional assay name.
#'
#' @return A data frame with one row per trial and M/H peak-to-peak
#'   amplitudes and fiducial samples.
#'
#' @references
#' Palmieri RM, Ingersoll CD, Hoffman MA (2004). The Hoffmann reflex:
#' methodologic considerations and applications for use in sports medicine and
#' athletic training research. *Journal of Athletic Training*, 39:268-277.
#'
#' @export
#'
#' @examples
#' pe <- make_mep(n_trials = 3, intensities = c(30, 40, 50), seed = 1)
#' hReflexRecruitment(
#'   pe, c(10, 20, 30),
#'   m_window_ms = c(10, 18), h_window_ms = c(20, 60),
#'   epoch_start_ms = -20
#' )
hReflexRecruitment <- function(x,
                               intensity,
                               m_window_ms = c(2, 12),
                               h_window_ms = c(20, 45),
                               epoch_start_ms = 0,
                               channel = 1,
                               assay_name = NULL) {
  source <- .mep_source(x, assay_name)
  epoch_start_ms <- .neuro_epoch_start(epoch_start_ms)
  if (!is.numeric(intensity) ||
      length(intensity) != source$n_trials ||
      any(!is.finite(intensity))) {
    stop("intensity must contain one finite numeric value per trial.",
         call. = FALSE)
  }
  selected <- .neuro_channel(source, channel)
  windows <- .neuro_response_windows(
    m_window_ms,
    h_window_ms,
    epoch_start_ms,
    source,
    "m_window_ms",
    "h_window_ms"
  )
  if (!length(selected) || source$n_trials == 0L) {
    return(.hreflex_empty())
  }

  rows <- vector("list", source$n_trials)
  for (trial in seq_len(source$n_trials)) {
    signal <- source$data[, selected, trial]
    m <- .mep_p2p(signal, windows$early)
    h <- .mep_p2p(signal, windows$late)
    rows[[trial]] <- data.frame(
      trial = as.integer(trial),
      intensity = as.numeric(intensity[trial]),
      m_p2p = unname(m["p2p"]),
      m_pos_peak_sample = as.integer(m["pos"]),
      m_neg_peak_sample = as.integer(m["neg"]),
      h_p2p = unname(h["p2p"]),
      h_pos_peak_sample = as.integer(h["pos"]),
      h_neg_peak_sample = as.integer(h["neg"])
    )
  }
  do.call(rbind, rows)
}


#' Calculate the Hmax/Mmax ratio
#'
#' @param recruitment An H-reflex recruitment table returned by
#'   [hReflexRecruitment()], or a data frame containing `intensity`, `h_p2p`,
#'   and `m_p2p`.
#'
#' @return A one-row data frame containing observed maxima, their lowest tied
#'   intensities, and the dimensionless Hmax/Mmax ratio.
#' @export
#'
#' @examples
#' recruitment <- data.frame(
#'   intensity = c(10, 20, 30),
#'   h_p2p = c(0.2, 1, 0.8),
#'   m_p2p = c(0.5, 1.5, 2)
#' )
#' hMaxMMax(recruitment)
hMaxMMax <- function(recruitment) {
  recruitment <- .hreflex_table(recruitment)
  if (!nrow(recruitment)) {
    return(data.frame(
      hmax = NA_real_,
      hmax_intensity = NA_real_,
      mmax = NA_real_,
      mmax_intensity = NA_real_,
      hmax_mmax_ratio = NA_real_
    ))
  }
  hmax <- max(recruitment$h_p2p)
  mmax <- max(recruitment$m_p2p)
  data.frame(
    hmax = hmax,
    hmax_intensity = min(recruitment$intensity[
      recruitment$h_p2p == hmax
    ]),
    mmax = mmax,
    mmax_intensity = min(recruitment$intensity[
      recruitment$m_p2p == mmax
    ]),
    hmax_mmax_ratio = if (mmax > 0) hmax / mmax else NA_real_
  )
}


#' Estimate H-reflex onset latency
#'
#' Uses the baseline-standard-deviation threshold and persistent-run detector
#' from [mepLatency()] within an H-reflex-specific search window.
#'
#' @inheritParams hReflexRecruitment
#' @param baseline_window_ms Optional baseline interval. The default uses all
#'   pre-stimulus samples.
#' @param threshold_sd Baseline standard-deviation multiplier.
#' @param min_consecutive_ms Required duration of a threshold crossing.
#'
#' @return A data frame containing onset sample and stimulus-relative H-reflex
#'   latency for each trial.
#' @export
#'
#' @examples
#' pe <- make_mep(n_trials = 3, mep_latency_ms = 25, seed = 2)
#' hReflexLatency(pe, epoch_start_ms = -20)
hReflexLatency <- function(x,
                           h_window_ms = c(20, 45),
                           epoch_start_ms = 0,
                           baseline_window_ms = NULL,
                           threshold_sd = 3,
                           min_consecutive_ms = 1,
                           channel = 1,
                           assay_name = NULL) {
  source <- .mep_source(x, assay_name)
  epoch_start_ms <- .neuro_epoch_start(epoch_start_ms)
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
  .mep_window(
    h_window_ms,
    epoch_start_ms,
    source$sr,
    source$n_time,
    "h_window_ms"
  )
  if (h_window_ms[1L] < 0) {
    stop("h_window_ms must be a post-stimulus window.", call. = FALSE)
  }
  if (!length(selected) || source$n_trials == 0L) {
    return(data.frame(
      trial = integer(),
      channel = integer(),
      channel_label = character(),
      onset_sample = integer(),
      h_latency_ms = numeric()
    ))
  }
  latency <- mepLatency(
    x,
    response_window_ms = h_window_ms,
    epoch_start_ms = epoch_start_ms,
    baseline_window_ms = baseline_window_ms,
    threshold_sd = threshold_sd,
    min_consecutive_ms = min_consecutive_ms,
    channels = selected,
    assay_name = assay_name
  )
  data.frame(
    trial = latency$trial,
    channel = latency$channel,
    channel_label = latency$channel_label,
    onset_sample = latency$onset_sample,
    h_latency_ms = latency$latency_ms
  )
}


#' Estimate the H-reflex recruitment threshold
#'
#' Finds the lowest stimulation intensity whose mean H amplitude reaches an
#' absolute cutoff or a fraction of the observed Mmax.
#'
#' @inheritParams hMaxMMax
#' @param criterion Non-negative absolute amplitude, or a fraction in `[0, 1]`
#'   when `scale = "mmax"`.
#' @param scale Whether `criterion` is relative to Mmax or in native amplitude
#'   units.
#'
#' @return A one-row data frame containing the first qualifying intensity,
#'   applied cutoff, and scale.
#' @export
#'
#' @examples
#' recruitment <- data.frame(
#'   intensity = c(10, 20, 30),
#'   h_p2p = c(0.1, 0.5, 1),
#'   m_p2p = c(0.2, 0.8, 2.5)
#' )
#' hReflexThreshold(recruitment, criterion = 0.2)
hReflexThreshold <- function(recruitment,
                             criterion = 0.05,
                             scale = c("mmax", "absolute")) {
  recruitment <- .hreflex_table(recruitment)
  scale <- match.arg(scale)
  upper <- if (scale == "mmax") 1 else Inf
  if (!is.numeric(criterion) || length(criterion) != 1L ||
      !is.finite(criterion) || criterion < 0 || criterion > upper) {
    stop(
      if (scale == "mmax") {
        "criterion must be one finite value in [0, 1] for scale = \"mmax\"."
      } else {
        "criterion must be one non-negative finite amplitude."
      },
      call. = FALSE
    )
  }

  if (!nrow(recruitment)) {
    return(data.frame(
      threshold_intensity = NA_real_,
      cutoff = if (scale == "absolute") criterion else NA_real_,
      scale = scale
    ))
  }
  mmax <- max(recruitment$m_p2p)
  cutoff <- if (scale == "mmax") criterion * mmax else criterion
  if (scale == "mmax" && mmax == 0) {
    threshold <- NA_real_
  } else {
    means <- stats::aggregate(
      recruitment$h_p2p,
      list(intensity = recruitment$intensity),
      mean
    )
    names(means)[2L] <- "h_mean"
    means <- means[order(means$intensity), , drop = FALSE]
    qualifying <- which(means$h_mean >= cutoff)
    threshold <- if (length(qualifying)) {
      means$intensity[qualifying[1L]]
    } else {
      NA_real_
    }
  }
  data.frame(
    threshold_intensity = threshold,
    cutoff = cutoff,
    scale = scale
  )
}
