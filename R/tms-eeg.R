# TMS-EEG pulse removal, source-informed cleaning, and TEP measurement

.tms_eeg_choice <- function(value, choices, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !value %in% choices) {
    stop(sprintf(
      "%s must be exactly one of: %s.",
      name, paste(shQuote(choices), collapse = ", ")
    ), call. = FALSE)
  }
  value
}


.tms_eeg_scalar <- function(value, name, lower = -Inf,
                            lower_open = FALSE, integer = FALSE) {
  valid <- is.numeric(value) && length(value) == 1L &&
    is.finite(value)
  if (valid && integer) {
    valid <- value >= -.Machine$integer.max &&
      value <= .Machine$integer.max &&
      value == floor(value)
  }
  if (valid) {
    valid <- if (lower_open) value > lower else value >= lower
  }
  if (!valid) {
    qualifier <- if (integer) "integer" else "numeric scalar"
    stop(sprintf("%s must be a finite %s in its documented range.",
                 name, qualifier), call. = FALSE)
  }
  if (integer) as.integer(value) else as.numeric(value)
}


.tms_eeg_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("%s must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}


.tms_eeg_source <- function(x, assay_name = NULL, epoch_start_ms = 0) {
  if (!inherits(x, "PhysioExperiment")) {
    stop("x must be a PhysioExperiment object.", call. = FALSE)
  }
  epoch_start_ms <- .tms_eeg_scalar(
    epoch_start_ms, "epoch_start_ms"
  )
  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }
  if (!is.character(assay_name) || length(assay_name) != 1L ||
      is.na(assay_name) || !nzchar(assay_name) ||
      !assay_name %in% SummarizedExperiment::assayNames(x)) {
    stop("assay_name must identify exactly one existing assay.",
         call. = FALSE)
  }

  original <- SummarizedExperiment::assay(x, assay_name)
  if (!is.numeric(original) || any(!is.finite(original))) {
    stop("The TMS-EEG assay must contain finite numeric data.",
         call. = FALSE)
  }
  dimensions <- dim(original)
  was_matrix <- length(dimensions) == 2L
  if (!was_matrix && length(dimensions) != 3L) {
    stop("The TMS-EEG assay must be time x channels or time x channels x trials.",
         call. = FALSE)
  }
  if (any(dimensions < 1L)) {
    stop("The TMS-EEG assay must contain samples, channels, and trials.",
         call. = FALSE)
  }

  data <- original
  if (was_matrix) {
    dn <- dimnames(original)
    if (is.null(dn)) {
      dn <- list(NULL, NULL)
    }
    data <- array(
      original,
      dim = c(dimensions, 1L),
      dimnames = c(dn, list(trial = "1"))
    )
  }
  dimensions3 <- dim(data)

  sr <- samplingRate(x)
  if (!is.numeric(sr) || length(sr) != 1L ||
      !is.finite(sr) || sr <= 0) {
    stop("The sampling rate must be a positive finite scalar.",
         call. = FALSE)
  }

  channel_data <- SummarizedExperiment::colData(x)
  if (!"label" %in% names(channel_data)) {
    stop("colData(x)$label is required for TMS-EEG analysis.",
         call. = FALSE)
  }
  labels <- as.character(channel_data$label)
  if (length(labels) != dimensions3[2L] || anyNA(labels) ||
      any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop("colData(x)$label must contain exact unique non-empty channel labels.",
         call. = FALSE)
  }

  list(
    data = data,
    original = original,
    original_dim = dim(original),
    original_dimnames = dimnames(original),
    was_matrix = was_matrix,
    assay_name = assay_name,
    sr = as.numeric(sr),
    epoch_start_ms = epoch_start_ms,
    n_time = dimensions3[1L],
    n_channels = dimensions3[2L],
    n_trials = dimensions3[3L],
    labels = labels
  )
}


.tms_eeg_window <- function(window_ms, source, name) {
  if (!is.numeric(window_ms) || length(window_ms) != 2L ||
      any(!is.finite(window_ms)) || window_ms[2L] <= window_ms[1L]) {
    stop(sprintf("%s must contain two increasing finite times.", name),
         call. = FALSE)
  }
  position <- round(
    (as.numeric(window_ms) - source$epoch_start_ms) / 1000 * source$sr
  ) + 1
  if (any(abs(position) > .Machine$integer.max)) {
    stop(sprintf("%s resolves outside the supported sample-index range.", name),
         call. = FALSE)
  }
  raw <- as.integer(position)
  if (raw[2L] < 1L || raw[1L] > source$n_time) {
    stop(sprintf("%s does not overlap the TMS-EEG epoch.", name),
         call. = FALSE)
  }
  range <- c(max(1L, raw[1L]), min(source$n_time, raw[2L]))
  if (range[2L] < range[1L]) {
    stop(sprintf("%s does not resolve to an epoch sample.", name),
         call. = FALSE)
  }
  list(
    window_ms = as.numeric(window_ms),
    raw_range = raw,
    range = range,
    samples = seq.int(range[1L], range[2L])
  )
}


.tms_eeg_windows <- function(pulse_window_ms, recharge_windows_ms, source) {
  if (!is.numeric(pulse_window_ms) || length(pulse_window_ms) != 2L ||
      any(!is.finite(pulse_window_ms)) ||
      pulse_window_ms[2L] <= pulse_window_ms[1L]) {
    stop("pulse_window_ms must contain two increasing finite numeric times.",
         call. = FALSE)
  }
  pulse <- matrix(as.numeric(pulse_window_ms), nrow = 1L)
  if (is.null(recharge_windows_ms)) {
    recharge <- matrix(numeric(), nrow = 0L, ncol = 2L)
  } else if (is.numeric(recharge_windows_ms) &&
             is.null(dim(recharge_windows_ms)) &&
             length(recharge_windows_ms) == 2L) {
    recharge <- matrix(as.numeric(recharge_windows_ms), nrow = 1L)
  } else if (is.numeric(recharge_windows_ms) &&
             length(dim(recharge_windows_ms)) == 2L &&
             ncol(recharge_windows_ms) == 2L) {
    recharge <- unname(as.matrix(recharge_windows_ms))
  } else {
    stop("recharge_windows_ms must be NULL, one pair, or an n x 2 matrix.",
         call. = FALSE)
  }
  requested <- rbind(pulse, recharge)
  if (any(!is.finite(requested)) ||
      any(requested[, 2L] <= requested[, 1L])) {
    stop("All pulse and recharge windows must be increasing and finite.",
         call. = FALSE)
  }
  resolved <- lapply(seq_len(nrow(requested)), function(i) {
    .tms_eeg_window(
      requested[i, ], source,
      if (i == 1L) "pulse_window_ms" else "recharge_windows_ms"
    )$range
  })
  ranges <- do.call(rbind, resolved)
  ranges <- ranges[order(ranges[, 1L], ranges[, 2L]), , drop = FALSE]
  merged <- list()
  for (i in seq_len(nrow(ranges))) {
    current <- as.integer(ranges[i, ])
    if (!length(merged) ||
        current[1L] > merged[[length(merged)]][2L] + 1L) {
      merged[[length(merged) + 1L]] <- current
    } else {
      merged[[length(merged)]][2L] <- max(
        merged[[length(merged)]][2L], current[2L]
      )
    }
  }
  merged <- do.call(rbind, merged)
  colnames(requested) <- c("start_ms", "end_ms")
  colnames(merged) <- c("first_sample", "last_sample")
  list(requested = requested, merged = merged)
}


.tms_eeg_restore <- function(data, source) {
  if (!source$was_matrix) {
    dimnames(data) <- dimnames(source$data)
    return(data)
  }
  matrix(
    data[, , 1L],
    nrow = source$n_time,
    ncol = source$n_channels,
    dimnames = source$original_dimnames
  )
}


.tms_eeg_flatten <- function(data) {
  dimensions <- dim(data)
  matrix(
    aperm(data, c(2L, 1L, 3L)),
    nrow = dimensions[2L],
    ncol = dimensions[1L] * dimensions[3L]
  )
}


.tms_eeg_write <- function(x, source, data, output_assay,
                           overwrite, step) {
  overwrite <- .tms_eeg_flag(overwrite, "overwrite")
  if (!is.character(output_assay) || length(output_assay) != 1L ||
      is.na(output_assay) || !nzchar(output_assay)) {
    stop("output_assay must be one non-empty name.", call. = FALSE)
  }
  existing <- SummarizedExperiment::assayNames(x)
  if (output_assay %in% existing &&
      !(identical(output_assay, source$assay_name) && overwrite)) {
    stop(
      "output_assay already exists; only the selected input may be explicitly overwritten.",
      call. = FALSE
    )
  }
  md <- S4Vectors::metadata(x)
  if (!is.null(md$tms_eeg) && !is.list(md$tms_eeg)) {
    stop("metadata(x)$tms_eeg must be a list.", call. = FALSE)
  }
  if (is.null(md$tms_eeg)) {
    md$tms_eeg <- list()
  }
  if (!is.null(md$tms_eeg$steps) && !is.list(md$tms_eeg$steps)) {
    stop("metadata(x)$tms_eeg$steps must be a list.", call. = FALSE)
  }
  if (is.null(md$tms_eeg$steps)) {
    md$tms_eeg$steps <- list()
  }

  restored <- .tms_eeg_restore(data, source)
  SummarizedExperiment::assay(x, output_assay) <- restored
  step$input_assay <- source$assay_name
  step$output_assay <- output_assay
  step$epoch_start_ms <- source$epoch_start_ms
  step$channels <- source$labels
  step$numeric_contract <- paste(
    "closed windows use round((time_ms-epoch_start_ms)/1000*sr)+1;",
    "time x channel x trial; exact channel identity"
  )
  md$tms_eeg$steps[[length(md$tms_eeg$steps) + 1L]] <- step
  S4Vectors::metadata(x) <- md
  x
}


.tms_eeg_mad <- function(x) {
  center <- stats::median(x)
  1.4826 * stats::median(abs(x - center))
}


.tms_eeg_slope <- function(values, samples) {
  centered <- samples - mean(samples)
  sum(centered * (values - mean(values))) / sum(centered^2)
}


.tms_eeg_pulse_provenance <- function(x, assay_name) {
  md <- S4Vectors::metadata(x)$tms_eeg
  steps <- md$steps
  if (is.list(steps) && length(steps)) {
    for (i in rev(seq_along(steps))) {
      step <- steps[[i]]
      if (is.list(step) && identical(step$method, "removeTMSpulse") &&
          identical(step$output_assay, assay_name)) {
        return(list(step = step, hardware_declared = FALSE))
      }
      if (is.list(step) && identical(step$output_assay, assay_name) &&
          is.list(step$pulse_provenance) &&
          (!is.null(step$pulse_provenance$step) ||
           identical(step$pulse_provenance$hardware_declared, TRUE))) {
        return(step$pulse_provenance)
      }
    }
  }
  list(
    step = NULL,
    hardware_declared = identical(md$pulse_absent, TRUE)
  )
}


.tms_eeg_pulse_step <- function(x, assay_name) {
  .tms_eeg_pulse_provenance(x, assay_name)$step
}


.tms_eeg_require_pulse_absent <- function(x, assay_name) {
  provenance <- .tms_eeg_pulse_provenance(x, assay_name)
  if (is.null(provenance$step) && !provenance$hardware_declared) {
    stop(
      "Remove the TMS pulse first or explicitly declare metadata(x)$tms_eeg$pulse_absent = TRUE.",
      call. = FALSE
    )
  }
  provenance
}


#' Remove TMS pulse and recharge intervals
#'
#' Replaces sample-exact closed pulse/recharge intervals by a linear or cubic
#' Hermite bridge, or blanks them to zero or a trial-specific baseline.
#' The selected input assay is preserved unless `overwrite = TRUE`.
#'
#' @details
#' Closed intervals use
#' `round((time_ms - epoch_start_ms) / 1000 * sampling_rate) + 1`.
#' Overlapping or sample-adjacent pulse/recharge intervals are merged. Linear
#' and cubic methods require real samples on both sides and never bridge across
#' a second removed interval. Baseline settings do not affect interpolation or
#' zero blanking. Baseline blanking requires at least two samples wholly outside
#' every removed interval. The audit metadata contains the complete altered
#' mask, merged sample ranges, continuity jumps, and robust baseline noise.
#'
#' @param x A `PhysioExperiment` with time x channel or time x channel x trial
#'   data.
#' @param pulse_window_ms Increasing pulse-window limits in milliseconds.
#' @param recharge_windows_ms Optional increasing pair or n x 2 matrix.
#' @param method Exactly `"cubic"`, `"linear"`, or `"blank"`.
#' @param support_ms One positive derivative-support duration for cubic
#'   interpolation, in milliseconds.
#' @param blank_value Exactly `"baseline"` or `"zero"`.
#' @param baseline_window_ms Closed baseline interval in milliseconds.
#' @param epoch_start_ms Time of the first sample relative to the pulse.
#' @param assay_name Optional input assay name.
#' @param output_assay Output assay name.
#' @param overwrite Whether to overwrite the selected input assay.
#'
#' @return A modified `PhysioExperiment` with an appended TMS-EEG audit step.
#'
#' @references
#' Rogasch NC, Sullivan C, Thomson RH, et al. (2017). Analysing concurrent
#' transcranial magnetic stimulation and electroencephalographic data.
#' *NeuroImage*, 147:934-951. \doi{10.1016/j.neuroimage.2016.10.031}
#'
#' @export
removeTMSpulse <- function(
    x,
    pulse_window_ms = c(-2, 10),
    recharge_windows_ms = NULL,
    method = c("cubic", "linear", "blank"),
    support_ms = 2,
    blank_value = c("baseline", "zero"),
    baseline_window_ms = c(-200, -20),
    epoch_start_ms = 0,
    assay_name = NULL,
    output_assay = "tms_pulse_clean",
    overwrite = FALSE) {
  if (missing(method)) {
    method <- method[1L]
  }
  if (length(method) != 1L) {
    stop("method must contain exactly one value.", call. = FALSE)
  }
  method <- .tms_eeg_choice(method, c("cubic", "linear", "blank"), "method")
  if (missing(blank_value)) {
    blank_value <- blank_value[1L]
  }
  if (length(blank_value) != 1L) {
    stop("blank_value must contain exactly one value.", call. = FALSE)
  }
  blank_value <- .tms_eeg_choice(
    blank_value, c("baseline", "zero"), "blank_value"
  )
  support_ms <- .tms_eeg_scalar(
    support_ms, "support_ms", lower = 0, lower_open = TRUE
  )
  source <- .tms_eeg_source(x, assay_name, epoch_start_ms)
  if (!is.numeric(baseline_window_ms) ||
      length(baseline_window_ms) != 2L ||
      any(!is.finite(baseline_window_ms)) ||
      baseline_window_ms[2L] <= baseline_window_ms[1L]) {
    stop("baseline_window_ms must contain two increasing finite times.",
         call. = FALSE)
  }
  windows <- .tms_eeg_windows(
    pulse_window_ms, recharge_windows_ms, source
  )
  ranges <- windows$merged
  if (any(ranges[, 1L] <= 1L) || any(ranges[, 2L] >= source$n_time)) {
    if (method != "blank") {
      stop("Every interpolated interval must leave a sample on both sides.",
           call. = FALSE)
    }
  }

  data <- source$data
  altered <- array(FALSE, dim(data), dimnames = dimnames(data))
  for (i in seq_len(nrow(ranges))) {
    altered[seq.int(ranges[i, 1L], ranges[i, 2L]), , ] <- TRUE
  }
  baseline_required <- method == "blank" && blank_value == "baseline"
  baseline_resolved <- if (baseline_required) {
    .tms_eeg_window(
      baseline_window_ms, source, "baseline_window_ms"
    )
  } else {
    tryCatch(
      .tms_eeg_window(
        baseline_window_ms, source, "baseline_window_ms"
      ),
      error = function(error) NULL
    )
  }
  baseline <- if (is.null(baseline_resolved)) {
    integer()
  } else {
    baseline_resolved$samples
  }
  baseline_clean <- baseline[!apply(
    altered[baseline, , , drop = FALSE], 1L, any
  )]
  if (baseline_required) {
    if (length(baseline_clean) != length(baseline) ||
        length(baseline_clean) < 2L) {
      stop(
        "Baseline blanking requires at least two baseline samples wholly outside removed intervals.",
        call. = FALSE
      )
    }
  }

  continuity <- vector("list", nrow(ranges) *
    source$n_channels * source$n_trials)
  row <- 0L
  support_samples <- as.integer(round(support_ms / 1000 * source$sr))
  for (gap in seq_len(nrow(ranges))) {
    first <- ranges[gap, 1L]
    last <- ranges[gap, 2L]
    left <- first - 1L
    right <- last + 1L
    previous_end <- if (gap == 1L) 0L else ranges[gap - 1L, 2L]
    next_start <- if (gap == nrow(ranges)) {
      source$n_time + 1L
    } else {
      ranges[gap + 1L, 1L]
    }
    if (method == "cubic") {
      left_support <- seq.int(
        max(previous_end + 1L, left - support_samples), left
      )
      right_support <- seq.int(
        right, min(next_start - 1L, right + support_samples)
      )
      if (length(left_support) < 2L || length(right_support) < 2L) {
        stop(
          "Cubic interpolation has insufficient real support beside an interval.",
          call. = FALSE
        )
      }
    }

    for (trial in seq_len(source$n_trials)) {
      for (channel in seq_len(source$n_channels)) {
        original <- source$data[, channel, trial]
        if (method == "linear") {
          samples <- seq.int(first, last)
          u <- (samples - left) / (right - left)
          data[samples, channel, trial] <-
            (1 - u) * original[left] + u * original[right]
        } else if (method == "cubic") {
          slope_left <- .tms_eeg_slope(
            original[left_support], left_support
          )
          slope_right <- .tms_eeg_slope(
            original[right_support], right_support
          )
          samples <- seq.int(first, last)
          width <- right - left
          u <- (samples - left) / width
          h00 <- 2 * u^3 - 3 * u^2 + 1
          h10 <- u^3 - 2 * u^2 + u
          h01 <- -2 * u^3 + 3 * u^2
          h11 <- u^3 - u^2
          data[samples, channel, trial] <-
            h00 * original[left] +
            h10 * width * slope_left +
            h01 * original[right] +
            h11 * width * slope_right
        } else {
          fill <- if (blank_value == "zero") {
            0
          } else {
            mean(original[baseline])
          }
          data[seq.int(first, last), channel, trial] <- fill
        }

        row <- row + 1L
        noise <- if (length(baseline_clean) >= 2L) {
          .tms_eeg_mad(original[baseline_clean])
        } else {
          NA_real_
        }
        continuity[[row]] <- data.frame(
          gap = gap,
          channel = channel,
          channel_label = source$labels[channel],
          trial = trial,
          left_jump = if (left >= 1L) {
            abs(data[first, channel, trial] - original[left])
          } else {
            NA_real_
          },
          right_jump = if (right <= source$n_time) {
            abs(original[right] - data[last, channel, trial])
          } else {
            NA_real_
          },
          baseline_robust_noise = noise
        )
      }
    }
  }
  continuity <- do.call(rbind, continuity)

  altered_out <- .tms_eeg_restore(altered, source)
  step <- list(
    method = "removeTMSpulse",
    windows_ms = list(
      pulse = as.numeric(pulse_window_ms),
      recharge = if (is.null(recharge_windows_ms)) {
        NULL
      } else {
        unname(windows$requested[-1L, , drop = FALSE])
      }
    ),
    sample_ranges = data.frame(
      first_sample = ranges[, 1L],
      last_sample = ranges[, 2L]
    ),
    parameters = list(
      replacement = method,
      support_ms = support_ms,
      support_ignored = method != "cubic",
      blank_value = blank_value,
      blank_settings_ignored = method != "blank",
      baseline_ignored = !baseline_required,
      baseline_window_ms = as.numeric(baseline_window_ms)
    ),
    diagnostics = list(
      altered_sample_mask = altered_out,
      continuity = continuity
    )
  )
  .tms_eeg_write(
    x, source, data, output_assay, overwrite, step
  )
}


.tms_eeg_rank <- function(x) {
  singular <- svd(x, nu = 0L, nv = 0L)$d
  if (!length(singular) || max(singular) == 0) {
    return(0L)
  }
  as.integer(sum(singular >
    max(dim(x)) * max(singular) * .Machine$double.eps))
}


.tms_eeg_leadfield <- function(forward_model, source, reference) {
  reference <- .tms_eeg_choice(reference, "average", "reference")
  if (!is.list(forward_model) || is.null(forward_model$leadfield)) {
    stop("forward_model must contain a numeric leadfield matrix.",
         call. = FALSE)
  }
  leadfield <- forward_model$leadfield
  if (!is.matrix(leadfield) || !is.numeric(leadfield) ||
      any(!is.finite(leadfield)) || ncol(leadfield) < 1L ||
      nrow(leadfield) != source$n_channels) {
    stop(
      "forward_model$leadfield must be finite with one row per EEG channel.",
      call. = FALSE
    )
  }

  row_labels <- rownames(leadfield)
  explicit <- forward_model$leadfield_channels
  if (is.null(explicit) && is.data.frame(forward_model$electrode_positions) &&
      "label" %in% names(forward_model$electrode_positions)) {
    explicit <- as.character(forward_model$electrode_positions$label)
  }
  valid_labels <- function(labels) {
    is.character(labels) && length(labels) == nrow(leadfield) &&
      !anyNA(labels) && all(nzchar(labels)) && !anyDuplicated(labels)
  }
  if (!is.null(explicit) && !valid_labels(explicit)) {
    stop("forward_model$leadfield_channels must contain exact unique labels.",
         call. = FALSE)
  }
  if (!is.null(row_labels) && valid_labels(row_labels) &&
      !is.null(explicit) &&
      !identical(as.character(row_labels), as.character(explicit))) {
    stop("Lead-field row names and leadfield_channels disagree.",
         call. = FALSE)
  }
  labels <- if (!is.null(explicit)) {
    as.character(explicit)
  } else {
    as.character(row_labels)
  }
  if (!valid_labels(labels) ||
      !setequal(labels, source$labels)) {
    stop("Lead-field labels must exactly cover the assay channel labels.",
         call. = FALSE)
  }
  order <- match(source$labels, labels)
  leadfield <- leadfield[order, , drop = FALSE]
  rownames(leadfield) <- source$labels

  channels <- source$n_channels
  H <- diag(channels) - matrix(1 / channels, channels, channels)
  referenced <- H %*% leadfield
  rank_h <- .tms_eeg_rank(H)
  rank_l <- .tms_eeg_rank(referenced)
  if (rank_l < 2L) {
    stop("The referenced lead field must have rank at least two.",
         call. = FALSE)
  }
  list(
    leadfield = leadfield,
    referenced = referenced,
    H = H,
    rank_H = rank_h,
    rank_leadfield = rank_l,
    labels = source$labels,
    method = forward_model$method,
    coordinate_provenance = forward_model$coordinate_provenance,
    reference = reference
  )
}


.tms_eeg_inverse <- function(matrix, name) {
  symmetric <- (matrix + t(matrix)) / 2
  eig <- eigen(symmetric, symmetric = TRUE)
  tolerance <- max(dim(symmetric)) * max(abs(eig$values)) *
    .Machine$double.eps
  if (min(eig$values) <= tolerance) {
    stop(sprintf("%s is singular or numerically indefinite.", name),
         call. = FALSE)
  }
  inverse <- sweep(eig$vectors, 2L, eig$values, "/") %*% t(eig$vectors)
  list(
    inverse = inverse,
    condition = max(eig$values) / min(eig$values),
    eigenvalues = eig$values
  )
}


.tms_eeg_orient <- function(u, labels) {
  for (component in seq_len(ncol(u))) {
    values <- abs(u[, component])
    candidates <- which(values == max(values))
    if (length(candidates) > 1L) {
      candidates <- candidates[order(
        enc2utf8(labels[candidates]), method = "radix"
      )]
    }
    if (u[candidates[1L], component] < 0) {
      u[, component] <- -u[, component]
    }
  }
  u
}


#' Clean TMS-EEG muscle artifact with SSP-SIR
#'
#' Estimates one deterministic artifact subspace across trials and reconstructs
#' the retained activity through a channel-matched lead field.
#'
#' @details
#' The pulse must already be removed or hardware absence must be explicitly
#' declared in `metadata(x)$tms_eeg$pulse_absent`. Both data and lead field are
#' average referenced with the same matrix. The returned data use the fixed
#' source-informed reconstruction `Lr G P`, not the plain SSP projection `P`.
#' SSP-SIR can attenuate neural signal when the artifact and neural subspaces
#' overlap; a lead field is an approximate physical model, not proof that a
#' retained deflection is cortical.
#'
#' @param x A pulse-cleaned `PhysioExperiment`.
#' @param forward_model A list containing a finite channel-matched `leadfield`.
#' @param artifact_window_ms Closed artifact-training interval.
#' @param n_artifact Optional exact number of artifact components.
#' @param variance_fraction Cumulative artifact variance target when
#'   `n_artifact` is `NULL`.
#' @param lambda Positive source-reconstruction regularization.
#' @param reference Currently exactly `"average"`.
#' @param epoch_start_ms Time of the first epoch sample relative to the pulse.
#' @param assay_name Optional pulse-cleaned input assay.
#' @param output_assay Output assay name.
#' @param overwrite Whether to overwrite the selected input assay.
#'
#' @return A modified `PhysioExperiment` with fixed SSP-SIR reconstruction
#'   metadata.
#'
#' @references
#' Mutanen TP, Kukkonen M, Nieminen JO, et al. (2016). Recovering TMS-evoked
#' EEG responses masked by muscle artifacts. *NeuroImage*, 139:157-166.
#' \doi{10.1016/j.neuroimage.2016.05.028}
#'
#' @export
tmsSSPSIR <- function(
    x,
    forward_model,
    artifact_window_ms = c(5, 50),
    n_artifact = NULL,
    variance_fraction = 0.9,
    lambda = 0.05,
    reference = "average",
    epoch_start_ms = 0,
    assay_name = NULL,
    output_assay = "tms_sspsir",
    overwrite = FALSE) {
  source <- .tms_eeg_source(x, assay_name, epoch_start_ms)
  pulse <- .tms_eeg_require_pulse_absent(x, source$assay_name)
  lead <- .tms_eeg_leadfield(forward_model, source, reference)
  lambda <- .tms_eeg_scalar(lambda, "lambda", 0, lower_open = TRUE)
  variance_fraction <- .tms_eeg_scalar(
    variance_fraction, "variance_fraction", 0, lower_open = TRUE
  )
  if (variance_fraction > 1) {
    stop("variance_fraction must be in (0, 1].", call. = FALSE)
  }
  if (!is.null(n_artifact)) {
    n_artifact <- .tms_eeg_scalar(
      n_artifact, "n_artifact", 1, integer = TRUE
    )
  }
  window <- .tms_eeg_window(
    artifact_window_ms, source, "artifact_window_ms"
  )
  if (window$range[1L] <= 1L) {
    stop("artifact_window_ms must leave a pre-window sample.",
         call. = FALSE)
  }

  flat <- .tms_eeg_flatten(source$data)
  referenced <- array(
    lead$H %*% flat,
    dim = c(source$n_channels, source$n_time, source$n_trials)
  )
  training <- lapply(seq_len(source$n_trials), function(trial) {
    trial_data <- referenced[, , trial, drop = FALSE][, , 1L]
    pre_mean <- rowMeans(
      trial_data[, seq_len(window$range[1L] - 1L), drop = FALSE]
    )
    sweep(
      trial_data[, window$samples, drop = FALSE],
      1L, pre_mean, "-"
    )
  })
  training <- do.call(cbind, training)
  decomposition <- svd(training, nu = source$n_channels, nv = 0L)
  u <- .tms_eeg_orient(decomposition$u, source$labels)
  fractions <- decomposition$d^2 / sum(decomposition$d^2)
  if (any(!is.finite(fractions))) {
    stop("The artifact training window has zero variance.",
         call. = FALSE)
  }
  selected <- if (is.null(n_artifact)) {
    which(cumsum(fractions) >= variance_fraction)[1L]
  } else {
    n_artifact
  }
  if (selected > ncol(u) ||
      source$n_channels - selected < 2L ||
      lead$rank_H - selected < 2L) {
    stop("Artifact components must leave at least two signal dimensions.",
         call. = FALSE)
  }

  Ua <- u[, seq_len(selected), drop = FALSE]
  P <- diag(source$n_channels) - Ua %*% t(Ua)
  Lr <- lead$referenced
  Lp <- P %*% Lr
  scale <- sum(Lp^2) / lead$rank_H
  if (!is.finite(scale) || scale <= 0) {
    stop("The projected lead field has no reconstructable energy.",
         call. = FALSE)
  }
  system <- Lp %*% t(Lp) +
    diag(lambda * scale, source$n_channels)
  solved <- .tms_eeg_inverse(system, "SSP-SIR regularized system")
  G <- t(Lp) %*% solved$inverse
  R_sir <- Lr %*% G %*% P
  clean_flat <- R_sir %*% lead$H %*% flat
  clean <- array(
    clean_flat,
    dim = c(source$n_channels, source$n_time, source$n_trials)
  )
  clean <- aperm(clean, c(2L, 1L, 3L))
  dimnames(clean) <- dimnames(source$data)
  plain_ssp <- P %*% lead$H %*% flat
  neural_error <- norm(R_sir %*% Lr - Lr, type = "F") /
    max(norm(Lr, type = "F"), .Machine$double.eps)

  step <- list(
    method = "tmsSSPSIR",
    windows_ms = list(artifact = as.numeric(artifact_window_ms)),
    sample_ranges = data.frame(
      first_sample = window$range[1L],
      last_sample = window$range[2L]
    ),
    parameters = list(
      n_artifact = selected,
      requested_n_artifact = n_artifact,
      variance_fraction = variance_fraction,
      lambda = lambda,
      reference = reference
    ),
    diagnostics = list(
      singular_values = decomposition$d,
      explained_fractions = fractions,
      artifact_topographies = Ua,
      projector = P,
      plain_ssp = plain_ssp,
      reconstruction_error = neural_error,
      effective_ranks = c(
        reference = lead$rank_H,
        leadfield = lead$rank_leadfield,
        projected_leadfield = .tms_eeg_rank(Lp)
      )
    ),
    lead_field = list(
      channels = lead$labels,
      dimensions = dim(lead$leadfield),
      method = lead$method,
      coordinate_provenance = lead$coordinate_provenance,
      reference_matrix = lead$H,
      condition_number = solved$condition,
      cleaning_matrix = R_sir
    ),
    pulse_provenance = pulse
  )
  .tms_eeg_write(
    x, source, clean, output_assay, overwrite, step
  )
}


#' Clean TMS-EEG sensor noise with SOUND
#'
#' Implements simultaneous leave-one-sensor noise estimation followed by the
#' source-space Wiener reconstruction of the SOUND algorithm.
#'
#' @details
#' SOUND assumes diagonal sensor noise, an approximate channel-matched lead
#' field, and i.i.d. source amplitudes for its regularized reconstruction. All
#' channel noise estimates are updated simultaneously from the preceding round,
#' and one final cleaning matrix is applied to every sample and trial. Pulse
#' samples are excluded from noise estimation. Regularization is not a
#' significance threshold and is not optimized against the reported TEP peaks.
#'
#' @inheritParams tmsSSPSIR
#' @param iterations Maximum number of complete simultaneous updates.
#' @param tolerance Non-negative maximum relative noise change for convergence.
#' @param analysis_window_ms Optional closed noise-estimation interval. `NULL`
#'   uses all samples except recorded pulse/recharge gaps.
#' @param output_assay Output assay name.
#'
#' @return A modified `PhysioExperiment` with the fixed SOUND cleaning matrix,
#'   complete noise history, and convergence diagnostics.
#'
#' @references
#' Mutanen TP, Metsomaa J, Liljander S, Ilmoniemi RJ (2018). Automatic and
#' robust noise suppression in EEG and MEG: The SOUND algorithm.
#' *NeuroImage*, 166:135-151. \doi{10.1016/j.neuroimage.2017.10.021}
#'
#' @export
soundClean <- function(
    x,
    forward_model,
    lambda = 0.1,
    iterations = 10L,
    tolerance = 0.01,
    reference = "average",
    analysis_window_ms = NULL,
    epoch_start_ms = 0,
    assay_name = NULL,
    output_assay = "tms_sound",
    overwrite = FALSE) {
  source <- .tms_eeg_source(x, assay_name, epoch_start_ms)
  if (source$n_channels < 3L) {
    stop("SOUND requires at least three channels.", call. = FALSE)
  }
  pulse <- .tms_eeg_require_pulse_absent(x, source$assay_name)
  lead <- .tms_eeg_leadfield(forward_model, source, reference)
  lambda <- .tms_eeg_scalar(lambda, "lambda", 0, lower_open = TRUE)
  iterations <- .tms_eeg_scalar(
    iterations, "iterations", 1, integer = TRUE
  )
  tolerance <- .tms_eeg_scalar(tolerance, "tolerance", 0)

  excluded <- integer()
  if (!is.null(pulse$step) && nrow(pulse$step$sample_ranges)) {
    excluded <- unique(unlist(lapply(
      seq_len(nrow(pulse$step$sample_ranges)),
      function(i) seq.int(
        pulse$step$sample_ranges$first_sample[i],
        pulse$step$sample_ranges$last_sample[i]
      )
    )))
  }
  if (is.null(analysis_window_ms)) {
    samples <- setdiff(seq_len(source$n_time), excluded)
    window_range <- NULL
  } else {
    window <- .tms_eeg_window(
      analysis_window_ms, source, "analysis_window_ms"
    )
    if (length(intersect(window$samples, excluded))) {
      stop("analysis_window_ms must exclude recorded pulse/recharge samples.",
           call. = FALSE)
    }
    samples <- window$samples
    window_range <- window$range
  }
  if (!length(samples)) {
    stop("No samples remain for SOUND noise estimation.", call. = FALSE)
  }

  flat <- .tms_eeg_flatten(source$data)
  referenced <- lead$H %*% flat
  selected_columns <- unlist(lapply(
    seq_len(source$n_trials),
    function(trial) samples + (trial - 1L) * source$n_time
  ))
  analysis <- referenced[, selected_columns, drop = FALSE]

  sigma <- rep(1, source$n_channels)
  history <- matrix(
    sigma, nrow = 1L,
    dimnames = list("initial", source$labels)
  )
  convergence <- numeric()
  diagnostics <- list()
  converged <- FALSE
  for (iteration in seq_len(iterations)) {
    next_sigma <- numeric(source$n_channels)
    for (sensor in seq_len(source$n_channels)) {
      minus <- setdiff(seq_len(source$n_channels), sensor)
      W_minus <- diag(1 / sigma[minus], length(minus))
      Lt <- W_minus %*% lead$referenced[minus, , drop = FALSE]
      yt <- W_minus %*% analysis[minus, , drop = FALSE]
      lambda_s <- lambda * sum(Lt^2) / length(minus)
      system <- Lt %*% t(Lt) + diag(lambda_s, length(minus))
      solved <- .tms_eeg_inverse(system, "SOUND leave-one-sensor system")
      J_minus <- t(Lt) %*% solved$inverse %*% yt
      predicted <- lead$referenced[sensor, , drop = FALSE] %*% J_minus
      residual <- as.numeric(
        analysis[sensor, ] - predicted
      )
      estimate <- sqrt(mean(residual^2))
      if (estimate == 0) {
        channel_rms <- sqrt(rowMeans(analysis^2))
        scale <- stats::median(channel_rms)
        estimate <- .Machine$double.eps *
          if (scale > 0) scale else 1
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(
          code = "zero_residual_clamped",
          iteration = iteration,
          channel = sensor,
          channel_label = source$labels[sensor]
        )
      }
      next_sigma[sensor] <- estimate
    }
    relative <- abs(next_sigma - sigma) / sigma
    convergence[iteration] <- max(relative)
    sigma <- next_sigma
    history <- rbind(history, sigma)
    rownames(history)[nrow(history)] <- paste0("iteration_", iteration)
    if (all(relative <= tolerance)) {
      converged <- TRUE
      break
    }
  }
  colnames(history) <- source$labels

  Wh <- diag(1 / sigma, source$n_channels)
  Lt <- Wh %*% lead$referenced
  lambda_final <- lambda * sum(Lt^2) / lead$rank_H
  system <- Lt %*% t(Lt) +
    diag(lambda_final, source$n_channels)
  solved <- .tms_eeg_inverse(system, "SOUND final Wiener system")
  K <- lead$referenced %*% t(Lt) %*%
    solved$inverse %*% Wh
  clean_flat <- K %*% referenced
  clean <- array(
    clean_flat,
    dim = c(source$n_channels, source$n_time, source$n_trials)
  )
  clean <- aperm(clean, c(2L, 1L, 3L))
  dimnames(clean) <- dimnames(source$data)

  typed_diagnostics <- if (length(diagnostics)) {
    do.call(rbind, diagnostics)
  } else {
    data.frame(
      code = character(),
      iteration = integer(),
      channel = integer(),
      channel_label = character()
    )
  }
  step <- list(
    method = "soundClean",
    windows_ms = list(analysis = analysis_window_ms),
    sample_ranges = if (is.null(window_range)) {
      data.frame(
        first_sample = min(samples),
        last_sample = max(samples)
      )
    } else {
      data.frame(
        first_sample = window_range[1L],
        last_sample = window_range[2L]
      )
    },
    parameters = list(
      lambda = lambda,
      iterations = iterations,
      tolerance = tolerance,
      reference = reference
    ),
    diagnostics = list(
      analysis_samples = samples,
      excluded_samples = excluded,
      noise_history = history,
      convergence_metric = convergence,
      converged = converged,
      completed_iterations = length(convergence),
      typed = typed_diagnostics
    ),
    lead_field = list(
      channels = lead$labels,
      dimensions = dim(lead$leadfield),
      method = lead$method,
      coordinate_provenance = lead$coordinate_provenance,
      reference_matrix = lead$H,
      rank = lead$rank_leadfield,
      condition_number = solved$condition,
      regularization = lambda_final,
      cleaning_matrix = K,
      final_noise_sd = sigma
    ),
    pulse_provenance = pulse
  )
  .tms_eeg_write(
    x, source, clean, output_assay, overwrite, step
  )
}


.tms_eeg_channels <- function(channels, source) {
  if (is.null(channels)) {
    return(seq_len(source$n_channels))
  }
  if (is.character(channels)) {
    if (anyNA(channels) || any(!nzchar(channels)) ||
        anyDuplicated(channels)) {
      stop("channels must contain exact unique labels.", call. = FALSE)
    }
    indices <- match(channels, source$labels)
    if (anyNA(indices)) {
      stop("Unknown TMS-EEG channel label.", call. = FALSE)
    }
    return(indices)
  }
  if (is.numeric(channels) && all(is.finite(channels)) &&
      all(channels >= -.Machine$integer.max) &&
      all(channels <= .Machine$integer.max) &&
      all(channels == floor(channels))) {
    indices <- as.integer(channels)
    if (!length(indices) || any(indices < 1L) ||
        any(indices > source$n_channels) || anyDuplicated(indices)) {
      stop("channels contain invalid or duplicate indices.",
           call. = FALSE)
    }
    return(indices)
  }
  stop("channels must contain exact labels or integer indices.",
       call. = FALSE)
}


.tms_eeg_component_windows <- function(components, windows_ms) {
  defaults <- rbind(
    N15 = c(10, 20),
    P30 = c(20, 40),
    N45 = c(35, 55),
    P60 = c(50, 75),
    N100 = c(80, 140),
    P180 = c(140, 220)
  )
  if (!is.character(components) || !length(components) ||
      anyNA(components) || anyDuplicated(components) ||
      any(!components %in% rownames(defaults))) {
    stop("components must be exact unique canonical TEP names.",
         call. = FALSE)
  }
  if (is.null(windows_ms)) {
    return(defaults[components, , drop = FALSE])
  }
  if (!is.numeric(windows_ms) ||
      length(dim(windows_ms)) != 2L || ncol(windows_ms) != 2L ||
      is.null(rownames(windows_ms)) ||
      anyNA(rownames(windows_ms)) || anyDuplicated(rownames(windows_ms)) ||
      !setequal(rownames(windows_ms), components)) {
    stop("windows_ms must be a named numeric matrix exactly covering components.",
         call. = FALSE)
  }
  windows <- windows_ms[components, , drop = FALSE]
  if (any(!is.finite(windows)) ||
      any(windows[, 2L] <= windows[, 1L])) {
    stop("Every component window must be increasing and finite.",
         call. = FALSE)
  }
  windows
}


#' Average TMS-evoked potentials and measure canonical peaks
#'
#' Baseline-corrects every trial and channel before a mean or trimmed mean,
#' then finds signed peaks independently in each requested channel.
#'
#' @details
#' Default signed peak windows are N15 (10--20 ms), P30 (20--40 ms), N45
#' (35--55 ms), P60 (50--75 ms), N100 (80--140 ms), and P180 (140--220 ms).
#' Exact ties use the earliest sample and a peak on a window endpoint is marked
#' as a boundary diagnostic. Canonical names define measurement windows only.
#' Auditory and somatosensory responses can remain, and no component is evidence
#' of a verified cortical source.
#'
#' @param x A `PhysioExperiment` containing epoched TMS-EEG.
#' @param components Exact canonical component names in requested order.
#' @param windows_ms Optional named component x 2 override matrix.
#' @param channels Optional exact labels or 1-based indices.
#' @param baseline_window_ms Closed baseline interval.
#' @param epoch_start_ms Time of the first sample relative to the pulse.
#' @param average Exactly `"mean"` or `"trimmed"`.
#' @param trim Trim fraction in `[0, 0.5)`.
#' @param assay_name Optional assay name.
#'
#' @return A `tms_tep` list with waveform, component, and settings tables.
#'
#' @references
#' Rogasch NC, Sullivan C, Thomson RH, et al. (2017). Analysing concurrent
#' transcranial magnetic stimulation and electroencephalographic data.
#' *NeuroImage*, 147:934-951. \doi{10.1016/j.neuroimage.2016.10.031}
#'
#' @export
tepAverage <- function(
    x,
    components = c("N15", "P30", "N45", "P60", "N100", "P180"),
    windows_ms = NULL,
    channels = NULL,
    baseline_window_ms = c(-200, -20),
    epoch_start_ms = 0,
    average = c("mean", "trimmed"),
    trim = 0.2,
    assay_name = NULL) {
  if (missing(average)) {
    average <- average[1L]
  }
  if (length(average) != 1L) {
    stop("average must contain exactly one value.", call. = FALSE)
  }
  average <- .tms_eeg_choice(average, c("mean", "trimmed"), "average")
  trim <- .tms_eeg_scalar(trim, "trim", 0)
  if (trim >= 0.5) {
    stop("trim must be in [0, 0.5).", call. = FALSE)
  }
  source <- .tms_eeg_source(x, assay_name, epoch_start_ms)
  selected <- .tms_eeg_channels(channels, source)
  windows <- .tms_eeg_component_windows(components, windows_ms)
  baseline <- .tms_eeg_window(
    baseline_window_ms, source, "baseline_window_ms"
  )
  if (length(baseline$samples) < 2L) {
    stop("baseline_window_ms must contain at least two samples.",
         call. = FALSE)
  }

  pulse <- .tms_eeg_pulse_step(x, source$assay_name)
  hardware <- identical(
    S4Vectors::metadata(x)$tms_eeg$pulse_absent, TRUE
  )
  if (is.null(pulse) && !hardware) {
    warning(
      "TMS pulse-cleaning provenance is missing; hardware blanking may have occurred.",
      call. = FALSE
    )
  }

  corrected <- source$data[, selected, , drop = FALSE]
  for (trial in seq_len(source$n_trials)) {
    trial_data <- matrix(
      corrected[, , trial, drop = FALSE],
      nrow = source$n_time,
      ncol = length(selected)
    )
    baseline_data <- trial_data[baseline$samples, , drop = FALSE]
    offsets <- colMeans(
      baseline_data
    )
    corrected[, , trial] <- sweep(
      trial_data,
      2L, offsets, "-"
    )
  }
  averaged <- apply(
    corrected, c(1L, 2L),
    function(values) {
      if (average == "mean") {
        mean(values)
      } else {
        mean(values, trim = trim)
      }
    }
  )
  if (length(selected) == 1L) {
    averaged <- matrix(averaged, ncol = 1L)
  }
  times <- source$epoch_start_ms +
    (seq_len(source$n_time) - 1L) / source$sr * 1000

  waveform <- do.call(rbind, lapply(seq_along(selected), function(i) {
    data.frame(
      time_ms = times,
      channel = selected[i],
      channel_label = source$labels[selected[i]],
      amplitude = averaged[, i],
      n_trials = source$n_trials
    )
  }))
  rownames(waveform) <- NULL

  polarity <- ifelse(substr(components, 1L, 1L) == "N",
                     "negative", "positive")
  component_rows <- vector("list", length(components) * length(selected))
  row <- 0L
  for (component in seq_along(components)) {
    resolved <- .tms_eeg_window(
      windows[component, ], source,
      paste0(components[component], " window")
    )
    for (i in seq_along(selected)) {
      values <- averaged[resolved$samples, i]
      local <- if (polarity[component] == "negative") {
        which.min(values)
      } else {
        which.max(values)
      }
      peak <- resolved$samples[local]
      row <- row + 1L
      component_rows[[row]] <- data.frame(
        component = components[component],
        polarity = polarity[component],
        window_start_ms = windows[component, 1L],
        window_end_ms = windows[component, 2L],
        channel = selected[i],
        channel_label = source$labels[selected[i]],
        amplitude = averaged[peak, i],
        latency_ms = times[peak],
        peak_sample = peak,
        boundary_peak = peak %in% resolved$range
      )
    }
  }
  component_table <- do.call(rbind, component_rows)
  rownames(component_table) <- NULL
  result <- list(
    waveform = waveform,
    components = component_table,
    settings = list(
      assay = source$assay_name,
      baseline_window_ms = as.numeric(baseline_window_ms),
      baseline_samples = baseline$range,
      average = average,
      trim = trim,
      pulse_cleaning = pulse,
      hardware_pulse_absent = hardware
    )
  )
  class(result) <- "tms_tep"
  result
}
