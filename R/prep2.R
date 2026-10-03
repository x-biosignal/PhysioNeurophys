# PREP2 upper-limb prognosis

#' Apply the PREP2 upper-limb prognosis algorithm
#'
#' Implements the day-3 PREP2 decision tree from SAFE score, age, motor-evoked
#' potential status, and NIH Stroke Scale score.
#'
#' @param safe_score Shoulder Abduction and Finger Extension MRC sum, an
#'   integer from 0 to 10.
#' @param age Age in years.
#' @param mep_status Optional MEP status: logical, `"present"`/`"MEP+"`, or
#'   `"absent"`/`"MEP-"`.
#' @param nihss Optional day-3 NIHSS total, an integer from 0 to 42.
#'
#' @return A `prep2_result` containing category, decision pathway, inputs, and
#'   a clinical description.
#'
#' @references
#' Stinear CM, Byblow WD, Ackerley SJ, Smith MC, Stinear JW, Barber PA
#' (2017). PREP2: A biomarker-based algorithm for predicting upper limb
#' function after stroke. *Annals of Clinical and Translational Neurology*,
#' 4:811-820. \doi{10.1002/acn3.488}
#'
#' @export
#'
#' @examples
#' prep2(safe_score = 8, age = 65)
#' prep2(2, 70, mep_status = "absent", nihss = 5)
prep2 <- function(safe_score, age, mep_status = NULL, nihss = NULL) {
  integer_range <- function(value, name, lower, upper) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value != as.integer(value) ||
        value < lower || value > upper) {
      stop(sprintf("%s must be one integer from %d to %d.",
                   name, lower, upper), call. = FALSE)
    }
    as.integer(value)
  }
  safe_score <- integer_range(safe_score, "safe_score", 0, 10)
  if (!is.numeric(age) || length(age) != 1L ||
      !is.finite(age) || age <= 0) {
    stop("age must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.null(nihss)) {
    nihss <- integer_range(nihss, "nihss", 0, 42)
  }
  normalize_mep <- function(value) {
    if (is.null(value)) {
      return(NULL)
    }
    if (is.logical(value) && length(value) == 1L && !is.na(value)) {
      return(if (value) "present" else "absent")
    }
    if (is.character(value) && length(value) == 1L && !is.na(value)) {
      key <- tolower(trimws(value))
      if (key %in% c("present", "mep+")) {
        return("present")
      }
      if (key %in% c("absent", "mep-")) {
        return("absent")
      }
    }
    stop(
      "mep_status must be present/MEP+, absent/MEP-, TRUE, or FALSE.",
      call. = FALSE
    )
  }
  mep_status <- normalize_mep(mep_status)

  if (safe_score >= 5L) {
    if (age < 80) {
      category <- "Excellent"
      pathway <- c("SAFE>=5", "age<80")
    } else {
      category <- "Good"
      pathway <- c("SAFE>=5", "age>=80")
    }
  } else {
    if (is.null(mep_status)) {
      stop("MEP status is required when SAFE < 5.", call. = FALSE)
    }
    if (mep_status == "present") {
      category <- "Good"
      pathway <- c("SAFE<5", "MEP+")
    } else {
      if (is.null(nihss)) {
        stop("Day-3 NIHSS is required when MEP is absent.",
             call. = FALSE)
      }
      if (nihss < 7L) {
        category <- "Limited"
        pathway <- c("SAFE<5", "MEP-", "NIHSS<7")
      } else {
        category <- "Poor"
        pathway <- c("SAFE<5", "MEP-", "NIHSS>=7")
      }
    }
  }
  descriptions <- c(
    Excellent = "Full or near-full upper-limb recovery is expected.",
    Good = "Substantial upper-limb recovery is expected, but not full recovery.",
    Limited = "Modest recovery is expected; emphasize task-specific training.",
    Poor = "Minimal upper-limb recovery is expected."
  )
  structure(
    list(
      category = category,
      pathway = pathway,
      inputs = list(
        safe_score = safe_score,
        age = as.numeric(age),
        mep_status = mep_status,
        nihss = nihss
      ),
      description = unname(descriptions[category])
    ),
    class = "prep2_result"
  )
}


#' @export
print.prep2_result <- function(x, ...) {
  cat("PREP2 prognosis:", x$category, "\n")
  cat("  Pathway:", paste(x$pathway, collapse = " -> "), "\n")
  cat(" ", x$description, "\n")
  invisible(x)
}
