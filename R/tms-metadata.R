# TMS pulse metadata

#' Empty TMS pulse table
#' @keywords internal
#' @noRd
.tms_empty <- function() {
  data.frame(
    pulse = integer(),
    intensity_pct_mso = numeric(),
    target = character(),
    coil = character(),
    rmt_pct_mso = numeric(),
    onset_sec = numeric()
  )
}


#' Resolve the trial count of a TMS experiment
#' @keywords internal
#' @noRd
.tms_trial_count <- function(x) {
  if (!inherits(x, "PhysioExperiment")) {
    stop("x must be a PhysioExperiment object.", call. = FALSE)
  }
  data <- SummarizedExperiment::assay(x, defaultAssay(x))
  dimensions <- dim(data)
  if (length(dimensions) == 2L) {
    1L
  } else if (length(dimensions) == 3L) {
    as.integer(dimensions[3L])
  } else {
    stop("The assay must be 2-D or epoched 3-D data.", call. = FALSE)
  }
}


#' Store TMS pulse metadata
#'
#' Stores one pulse row per trial without replacing other experiment metadata.
#' Intensity and resting motor threshold are expressed as percent of maximum
#' stimulator output (`%MSO`).
#'
#' @param x A `PhysioExperiment` with a 2-D or epoched 3-D assay.
#' @param intensity_pct_mso Numeric stimulation intensity, one value per trial.
#' @param target Cortical target or target-muscle label.
#' @param coil Coil descriptor.
#' @param rmt_pct_mso Resting motor threshold in `%MSO`.
#' @param onset_sec Pulse onset in the source continuous recording, in seconds.
#'
#' @return The modified `PhysioExperiment`.
#' @export
#'
#' @examples
#' pe <- make_mep(n_trials = 3, intensities = c(40, 50, 60))
#' pe <- setTMSpulses(pe, c(40, 50, 60), target = "M1_FDI")
#' getTMSpulses(pe)
setTMSpulses <- function(x,
                         intensity_pct_mso,
                         target = NA_character_,
                         coil = NA_character_,
                         rmt_pct_mso = NA_real_,
                         onset_sec = NA_real_) {
  n_trials <- .tms_trial_count(x)
  if (!is.numeric(intensity_pct_mso) ||
      length(intensity_pct_mso) != n_trials ||
      any(!is.finite(intensity_pct_mso)) ||
      any(intensity_pct_mso < 0 | intensity_pct_mso > 100)) {
    stop(
      "intensity_pct_mso must contain one finite value in [0, 100] per trial.",
      call. = FALSE
    )
  }

  recycle_character <- function(value, name) {
    if (!is.character(value) || !length(value) ||
        !(length(value) %in% c(1L, n_trials))) {
      stop(sprintf("%s must be character and scalar or one value per trial.",
                   name), call. = FALSE)
    }
    rep(value, length.out = n_trials)
  }
  recycle_numeric <- function(value, name, lower = -Inf, upper = Inf) {
    if (!is.numeric(value) || !length(value) ||
        !(length(value) %in% c(1L, n_trials)) ||
        any(!is.na(value) &
              (!is.finite(value) | value < lower | value > upper))) {
      stop(sprintf(
        "%s must be numeric, scalar or one value per trial, and within range.",
        name
      ), call. = FALSE)
    }
    rep(as.numeric(value), length.out = n_trials)
  }

  pulses <- data.frame(
    pulse = seq_len(n_trials),
    intensity_pct_mso = as.numeric(intensity_pct_mso),
    target = recycle_character(target, "target"),
    coil = recycle_character(coil, "coil"),
    rmt_pct_mso = recycle_numeric(
      rmt_pct_mso, "rmt_pct_mso", lower = 0, upper = 100
    ),
    onset_sec = recycle_numeric(onset_sec, "onset_sec", lower = 0)
  )
  metadata <- S4Vectors::metadata(x)
  metadata$tms <- pulses
  S4Vectors::metadata(x) <- metadata
  x
}


#' Retrieve TMS pulse metadata
#'
#' @param x A `PhysioExperiment`.
#'
#' @return A six-column data frame with one row per stored pulse, or a typed
#'   zero-row table when TMS metadata are absent.
#' @export
#'
#' @examples
#' getTMSpulses(make_mep(n_trials = 3))
getTMSpulses <- function(x) {
  if (!inherits(x, "PhysioExperiment")) {
    stop("x must be a PhysioExperiment object.", call. = FALSE)
  }
  pulses <- S4Vectors::metadata(x)$tms
  if (is.null(pulses)) {
    return(.tms_empty())
  }
  required <- names(.tms_empty())
  if (!is.data.frame(pulses) || !identical(names(pulses), required)) {
    stop("Stored TMS metadata do not follow the six-column schema.",
         call. = FALSE)
  }
  pulses
}
