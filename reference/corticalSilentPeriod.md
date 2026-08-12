# Measure the cortical silent period

Detects sustained suppression and recovery of rectified tonic EMG after
a TMS pulse using a baseline-relative threshold.

## Usage

``` r
corticalSilentPeriod(
  x,
  search_window_ms = c(50, 400),
  epoch_start_ms = 0,
  baseline_window_ms = NULL,
  threshold_fraction = 0.25,
  min_silence_ms = 10,
  min_return_ms = 5,
  channel = 1,
  assay_name = NULL
)
```

## Arguments

- x:

  An epoched tonic-EMG `PhysioExperiment`.

- search_window_ms:

  Post-stimulus CSP search window.

- epoch_start_ms:

  Time of the first epoch sample relative to the pulse.

- baseline_window_ms:

  Optional quiet-reference interval. The default uses all pre-stimulus
  samples.

- threshold_fraction:

  Fraction of median rectified baseline below which EMG is considered
  silent.

- min_silence_ms:

  Required duration of suppression.

- min_return_ms:

  Required duration of recovered activity.

- channel:

  One channel label or 1-based index.

- assay_name:

  Optional assay name.

## Value

A data frame with onset, offset, and duration for each trial.

## References

Rossini PM, Burke D, Chen R, et al. (2015). Non-invasive electrical and
magnetic stimulation of the brain, spinal cord, roots and peripheral
nerves: basic principles and procedures for routine clinical and
research application. *Clinical Neurophysiology*, 126:1071-1107.
[doi:10.1016/j.clinph.2015.02.001](https://doi.org/10.1016/j.clinph.2015.02.001)

## Examples

``` r
data <- array(rep(c(-1, 1), 250), c(500, 1, 1))
data[181:330, 1, 1] <- 0
pe <- PhysioExperiment(
  assays = list(raw = data),
  samplingRate = 1000
)
corticalSilentPeriod(pe, epoch_start_ms = -100)
#>   trial channel channel_label baseline_rectified onset_sample offset_sample
#> 1     1       1           Ch1                  1          181           331
#>   onset_ms offset_ms duration_ms
#> 1       80       230         150
```
