# Motor-evoked-potential analysis

#' Resolve and validate an MEP assay
#' @keywords internal
#' @noRd
.mep_source <- function(x, assay_name = NULL) {
  if (!inherits(x, "PhysioExperiment")) {
    stop("x must be a PhysioExperiment object.", call. = FALSE)
  }
  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }
  data <- SummarizedExperiment::assay(x, assay_name)
  if (!is.numeric(data) || any(!is.finite(data))) {
    stop("The MEP assay must contain finite numeric data.", call. = FALSE)
  }
  dimensions <- dim(data)
  if (length(dimensions) == 2L) {
    data_dimnames <- dimnames(data)
    if (is.null(data_dimnames)) {
      data_dimnames <- list(NULL, NULL)
    }
    data <- array(
      data,
      dim = c(dimensions[1L], dimensions[2L], 1L),
      dimnames = c(data_dimnames, list(trial = "1"))
    )
  } else if (length(dimensions) != 3L) {
    stop("The MEP assay must be time x channels or time x channels x trials.",
         call. = FALSE)
  }
  dimensions <- dim(data)
  if (dimensions[1L] < 1L) {
    stop("The MEP assay must contain at least one time sample.",
         call. = FALSE)
  }
  sr <- samplingRate(x)
  if (!is.numeric(sr) || length(sr) != 1L ||
      !is.finite(sr) || sr <= 0) {
    stop("The sampling rate must be a positive finite scalar.",
         call. = FALSE)
  }
  channel_data <- SummarizedExperiment::colData(x)
  labels <- if ("label" %in% names(channel_data)) {
    as.character(channel_data$label)
  } else if (dimensions[2L] == 0L) {
    character()
  } else {
    paste0("Ch", seq_len(dimensions[2L]))
  }
  if (length(labels) != dimensions[2L]) {
    stop("colData must contain one row per MEP channel.", call. = FALSE)
  }
  list(
    data = data,
    sr = as.numeric(sr),
    n_time = dimensions[1L],
    n_channels = dimensions[2L],
    n_trials = dimensions[3L],
    labels = labels
  )
}


#' Resolve channel selections
#' @keywords internal
#' @noRd
.mep_channels <- function(channels, labels, n_channels) {
  if (is.null(channels)) {
    return(seq_len(n_channels))
  }
  if (is.character(channels)) {
    indices <- match(channels, labels)
    if (anyNA(indices)) {
      stop("Unknown MEP channel label.", call. = FALSE)
    }
  } else if (is.numeric(channels) &&
             all(is.finite(channels)) &&
             all(channels == as.integer(channels))) {
    indices <- as.integer(channels)
  } else {
    stop("channels must contain channel labels or integer indices.",
         call. = FALSE)
  }
  if (any(indices < 1L | indices > n_channels) ||
      anyDuplicated(indices)) {
    stop("channels contain invalid or duplicate indices.", call. = FALSE)
  }
  indices
}


#' Convert a relative-time window to assay samples
#' @keywords internal
#' @noRd
.mep_window <- function(window_ms, epoch_start_ms, sr, n_time, name) {
  if (!is.numeric(window_ms) || length(window_ms) != 2L ||
      any(!is.finite(window_ms)) || window_ms[2L] <= window_ms[1L]) {
    stop(sprintf("%s must contain two increasing finite times.", name),
         call. = FALSE)
  }
  lo <- as.integer(round(
    (window_ms[1L] - epoch_start_ms) / 1000 * sr
  )) + 1L
  hi <- as.integer(round(
    (window_ms[2L] - epoch_start_ms) / 1000 * sr
  )) + 1L
  lo <- max(1L, lo)
  hi <- min(n_time, hi)
  if (hi < lo) {
    stop(sprintf("%s does not overlap the MEP epoch.", name),
         call. = FALSE)
  }
  seq.int(lo, hi)
}


#' Resolve per-trial stimulation intensity
#' @keywords internal
#' @noRd
.mep_intensity <- function(x, n_trials) {
  pulses <- getTMSpulses(x)
  if (nrow(pulses) == n_trials) {
    pulses$intensity_pct_mso
  } else {
    rep(NA_real_, n_trials)
  }
}


#' Peak-to-peak MEP fiducials
#' @keywords internal
#' @noRd
.mep_p2p <- function(signal, samples) {
  segment <- signal[samples]
  pos <- samples[which.max(segment)]
  neg <- samples[which.min(segment)]
  c(
    p2p = signal[pos] - signal[neg],
    pos = pos,
    neg = neg
  )
}


#' Empty MEP amplitude table
#' @keywords internal
#' @noRd
.mep_empty_amplitude <- function() {
  data.frame(
    channel = integer(),
    channel_label = character(),
    trial = integer(),
    intensity_pct_mso = numeric(),
    p2p = numeric(),
    pos_peak_sample = integer(),
    neg_peak_sample = integer(),
    pos_peak_time_ms = numeric(),
    neg_peak_time_ms = numeric()
  )
}


#' Measure motor-evoked-potential amplitude
#'
#' Calculates peak-to-peak MEP amplitude within a response window for every
#' selected channel and trial.
#'
#' @param x A `PhysioExperiment` with time x channels x trials MEP data.
#' @param response_window_ms Two response-window limits relative to the TMS
#'   pulse, in milliseconds.
#' @param epoch_start_ms Time of the first epoch sample relative to the pulse.
#' @param channels Optional channel labels or 1-based indices.
#' @param assay_name Optional assay name.
#'
#' @return A data frame with one row per selected channel and trial.
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
#' pe <- make_mep(n_trials = 4, seed = 1)
#' mepAmplitude(pe, epoch_start_ms = -20)
mepAmplitude <- function(x,
                         response_window_ms = c(10, 60),
                         epoch_start_ms = 0,
                         channels = NULL,
                         assay_name = NULL) {
  source <- .mep_source(x, assay_name)
  if (!is.numeric(epoch_start_ms) || length(epoch_start_ms) != 1L ||
      !is.finite(epoch_start_ms)) {
    stop("epoch_start_ms must be a finite numeric scalar.", call. = FALSE)
  }
  selected <- .mep_channels(
    channels, source$labels, source$n_channels
  )
  if (!length(selected) || source$n_trials == 0L) {
    return(.mep_empty_amplitude())
  }
  samples <- .mep_window(
    response_window_ms,
    epoch_start_ms,
    source$sr,
    source$n_time,
    "response_window_ms"
  )
  intensity <- .mep_intensity(x, source$n_trials)
  rows <- vector("list", length(selected) * source$n_trials)
  row <- 0L
  for (channel in selected) {
    for (trial in seq_len(source$n_trials)) {
      row <- row + 1L
      fiducials <- .mep_p2p(
        source$data[, channel, trial], samples
      )
      pos <- as.integer(fiducials["pos"])
      neg <- as.integer(fiducials["neg"])
      rows[[row]] <- data.frame(
        channel = as.integer(channel),
        channel_label = source$labels[channel],
        trial = as.integer(trial),
        intensity_pct_mso = intensity[trial],
        p2p = unname(fiducials["p2p"]),
        pos_peak_sample = pos,
        neg_peak_sample = neg,
        pos_peak_time_ms = epoch_start_ms +
          (pos - 1) / source$sr * 1000,
        neg_peak_time_ms = epoch_start_ms +
          (neg - 1) / source$sr * 1000
      )
    }
  }
  do.call(rbind, rows)
}


#' First persistent threshold crossing
#' @keywords internal
#' @noRd
.mep_onset <- function(exceed, samples, min_consecutive) {
  if (!length(exceed)) {
    return(NA_integer_)
  }
  runs <- rle(exceed)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  accepted <- which(runs$values & runs$lengths >= min_consecutive)
  if (!length(accepted)) {
    NA_integer_
  } else {
    as.integer(samples[starts[accepted[1L]]])
  }
}


#' Empty MEP latency table
#' @keywords internal
#' @noRd
.mep_empty_latency <- function() {
  data.frame(
    channel = integer(),
    channel_label = character(),
    trial = integer(),
    intensity_pct_mso = numeric(),
    onset_sample = integer(),
    latency_ms = numeric()
  )
}


#' Estimate motor-evoked-potential onset latency
#'
#' Detects the first persistent signed deviation from the baseline mean above
#' a multiple of the baseline standard deviation.
#'
#' @inheritParams mepAmplitude
#' @param response_window_ms Search window relative to the TMS pulse.
#' @param baseline_window_ms Optional baseline interval. The default uses all
#'   pre-stimulus samples.
#' @param threshold_sd Baseline standard-deviation multiplier.
#' @param min_consecutive_ms Required duration of a threshold crossing.
#'
#' @return A data frame containing onset sample and latency for each selected
#'   channel and trial.
#' @export
#'
#' @examples
#' pe <- make_mep(n_trials = 4, intensities = rep(65, 4), seed = 2)
#' mepLatency(pe, epoch_start_ms = -20)
mepLatency <- function(x,
                       response_window_ms = c(5, 50),
                       epoch_start_ms = 0,
                       baseline_window_ms = NULL,
                       threshold_sd = 3,
                       min_consecutive_ms = 1,
                       channels = NULL,
                       assay_name = NULL) {
  source <- .mep_source(x, assay_name)
  if (!is.numeric(epoch_start_ms) || length(epoch_start_ms) != 1L ||
      !is.finite(epoch_start_ms)) {
    stop("epoch_start_ms must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(threshold_sd) || length(threshold_sd) != 1L ||
      !is.finite(threshold_sd) || threshold_sd <= 0 ||
      !is.numeric(min_consecutive_ms) ||
      length(min_consecutive_ms) != 1L ||
      !is.finite(min_consecutive_ms) || min_consecutive_ms <= 0) {
    stop("threshold_sd and min_consecutive_ms must be positive scalars.",
         call. = FALSE)
  }
  selected <- .mep_channels(
    channels, source$labels, source$n_channels
  )
  if (!length(selected) || source$n_trials == 0L) {
    return(.mep_empty_latency())
  }
  response <- .mep_window(
    response_window_ms,
    epoch_start_ms,
    source$sr,
    source$n_time,
    "response_window_ms"
  )
  times <- epoch_start_ms +
    (seq_len(source$n_time) - 1) / source$sr * 1000
  warned_no_pre <- FALSE
  if (is.null(baseline_window_ms)) {
    baseline <- which(times < 0)
    if (length(baseline) < 2L) {
      baseline <- seq_len(min(
        source$n_time,
        max(2L, as.integer(round(0.020 * source$sr)))
      ))
      warning(
        "No pre-stimulus baseline; using the first 20 ms of the epoch.",
        call. = FALSE
      )
      warned_no_pre <- TRUE
    }
  } else {
    if (!is.numeric(baseline_window_ms) ||
        length(baseline_window_ms) != 2L ||
        any(!is.finite(baseline_window_ms)) ||
        baseline_window_ms[2L] <= baseline_window_ms[1L]) {
      stop("baseline_window_ms must contain two increasing finite times.",
           call. = FALSE)
    }
    baseline <- which(
      times >= baseline_window_ms[1L] &
        times < baseline_window_ms[2L]
    )
  }
  stim_sample <- as.integer(round(
    (0 - epoch_start_ms) / 1000 * source$sr
  )) + 1L
  search <- response[response >= stim_sample]
  if (!length(search)) {
    stop("response_window_ms contains no post-stimulus samples.",
         call. = FALSE)
  }
  min_consecutive <- max(
    1L,
    as.integer(round(min_consecutive_ms / 1000 * source$sr))
  )
  intensity <- .mep_intensity(x, source$n_trials)
  rows <- vector("list", length(selected) * source$n_trials)
  warned_degenerate <- FALSE
  row <- 0L
  for (channel in selected) {
    for (trial in seq_len(source$n_trials)) {
      row <- row + 1L
      signal <- source$data[, channel, trial]
      baseline_mean <- if (length(baseline)) {
        mean(signal[baseline])
      } else {
        mean(signal[search])
      }
      baseline_sd <- if (length(baseline) >= 2L) {
        stats::sd(signal[baseline])
      } else {
        NA_real_
      }
      if (!is.finite(baseline_sd) || baseline_sd == 0) {
        baseline_sd <- stats::sd(signal[search])
        if (!warned_degenerate && !warned_no_pre) {
          warning(
            "Degenerate baseline; using the response-window standard deviation.",
            call. = FALSE
          )
          warned_degenerate <- TRUE
        }
      }
      if (!is.finite(baseline_sd)) {
        baseline_sd <- 0
      }
      exceed <- abs(signal[search] - baseline_mean) >
        threshold_sd * baseline_sd
      onset <- .mep_onset(exceed, search, min_consecutive)
      rows[[row]] <- data.frame(
        channel = as.integer(channel),
        channel_label = source$labels[channel],
        trial = as.integer(trial),
        intensity_pct_mso = intensity[trial],
        onset_sample = onset,
        latency_ms = if (is.na(onset)) {
          NA_real_
        } else {
          epoch_start_ms + (onset - 1) / source$sr * 1000
        }
      )
    }
  }
  do.call(rbind, rows)
}


#' Fit a Boltzmann MEP recruitment curve
#'
#' Fits peak-to-peak MEP amplitude against stimulator intensity with a
#' three-parameter logistic curve.
#'
#' @param x An epoched MEP `PhysioExperiment`.
#' @param channel One channel label or 1-based index.
#' @param response_window_ms MEP response window in milliseconds.
#' @param epoch_start_ms Time of the first epoch sample relative to the pulse.
#' @param intensity Optional per-trial stimulation intensity in `%MSO`.
#' @param assay_name Optional assay name.
#'
#' @return A `recruitment_curve` object with fitted parameters, predictions,
#'   fit quality, and the underlying `nls` model.
#'
#' @references
#' Devanne H, Lavoie BA, Capaday C (1997). Input-output properties and gain
#' changes in the human corticospinal pathway. *Experimental Brain Research*,
#' 114:329-338. \doi{10.1007/PL00005641}
#'
#' @export
#'
#' @examples
#' pe <- make_mep(n_trials = 15, baseline_sd = 0.001, seed = 3)
#' recruitmentCurve(pe, epoch_start_ms = -20)
recruitmentCurve <- function(x,
                             channel = 1,
                             response_window_ms = c(10, 60),
                             epoch_start_ms = 0,
                             intensity = NULL,
                             assay_name = NULL) {
  source <- .mep_source(x, assay_name)
  selected <- .mep_channels(
    channel, source$labels, source$n_channels
  )
  if (length(selected) != 1L) {
    stop("channel must select exactly one MEP channel.", call. = FALSE)
  }
  amplitudes <- mepAmplitude(
    x,
    response_window_ms = response_window_ms,
    epoch_start_ms = epoch_start_ms,
    channels = selected,
    assay_name = assay_name
  )
  if (is.null(intensity)) {
    intensity <- amplitudes$intensity_pct_mso
  }
  if (!is.numeric(intensity) ||
      length(intensity) != source$n_trials ||
      any(!is.finite(intensity))) {
    stop("intensity must contain one finite numeric value per trial.",
         call. = FALSE)
  }
  if (length(unique(intensity)) < 3L) {
    stop("At least three distinct stimulation intensities are required.",
         call. = FALSE)
  }
  data <- data.frame(
    intensity = as.numeric(intensity),
    p2p = amplitudes$p2p
  )
  data <- data[order(data$intensity), , drop = FALSE]

  model <- tryCatch(
    suppressWarnings(stats::nls(
      p2p ~ stats::SSlogis(intensity, Asym, xmid, scal),
      data = data
    )),
    error = function(e) NULL
  )
  if (is.null(model)) {
    starts <- list(
      Asym = max(data$p2p),
      xmid = stats::median(data$intensity),
      scal = (max(data$intensity) - min(data$intensity)) / 4
    )
    model <- tryCatch(
      suppressWarnings(stats::nls(
        p2p ~ stats::SSlogis(intensity, Asym, xmid, scal),
        data = data,
        start = starts,
        control = stats::nls.control(maxiter = 200, warnOnly = TRUE)
      )),
      error = function(e) NULL
    )
  }

  if (is.null(model)) {
    plateau <- s50 <- k <- NA_real_
    data$fitted <- NA_real_
    converged <- FALSE
  } else {
    coefficients <- stats::coef(model)
    plateau <- unname(coefficients["Asym"])
    s50 <- unname(coefficients["xmid"])
    k <- unname(coefficients["scal"])
    data$fitted <- as.numeric(stats::predict(model, newdata = data))
    converged <- isTRUE(model$convInfo$isConv)
  }
  total <- sum((data$p2p - mean(data$p2p))^2)
  r_squared <- if (!is.finite(total) || total == 0 ||
                   any(!is.finite(data$fitted))) {
    NA_real_
  } else {
    1 - sum((data$p2p - data$fitted)^2) / total
  }
  max_slope <- if (is.finite(plateau) && is.finite(k) && k != 0) {
    plateau / (4 * k)
  } else {
    NA_real_
  }

  structure(
    list(
      plateau = plateau,
      s50 = s50,
      k = k,
      max_slope = max_slope,
      r_squared = r_squared,
      data = data,
      n = as.integer(nrow(data)),
      model = model,
      converged = converged
    ),
    class = "recruitment_curve"
  )
}


#' @export
print.recruitment_curve <- function(x, ...) {
  cat("MEP recruitment curve\n")
  cat(sprintf("  MEP_max: %s\n", format(x$plateau, digits = 5)))
  cat(sprintf("  S50: %s %%MSO\n", format(x$s50, digits = 5)))
  cat(sprintf("  k: %s %%MSO\n", format(x$k, digits = 5)))
  cat(sprintf("  Maximum slope: %s\n", format(x$max_slope, digits = 5)))
  cat(sprintf("  R-squared: %s\n", format(x$r_squared, digits = 5)))
  cat(sprintf("  Trials: %d; converged: %s\n", x$n, x$converged))
  invisible(x)
}
