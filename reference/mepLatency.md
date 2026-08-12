# Estimate motor-evoked-potential onset latency

Detects the first persistent signed deviation from the baseline mean
above a multiple of the baseline standard deviation.

## Usage

``` r
mepLatency(
  x,
  response_window_ms = c(5, 50),
  epoch_start_ms = 0,
  baseline_window_ms = NULL,
  threshold_sd = 3,
  min_consecutive_ms = 1,
  channels = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A `PhysioExperiment` with time x channels x trials MEP data.

- response_window_ms:

  Search window relative to the TMS pulse.

- epoch_start_ms:

  Time of the first epoch sample relative to the pulse.

- baseline_window_ms:

  Optional baseline interval. The default uses all pre-stimulus samples.

- threshold_sd:

  Baseline standard-deviation multiplier.

- min_consecutive_ms:

  Required duration of a threshold crossing.

- channels:

  Optional channel labels or 1-based indices.

- assay_name:

  Optional assay name.

## Value

A data frame containing onset sample and latency for each selected
channel and trial.

## Examples

``` r
pe <- make_mep(n_trials = 4, intensities = rep(65, 4), seed = 2)
mepLatency(pe, epoch_start_ms = -20)
#>   channel channel_label trial intensity_pct_mso onset_sample latency_ms
#> 1       1          MEP1     1                65          202       20.2
#> 2       1          MEP1     2                65          202       20.2
#> 3       1          MEP1     3                65          202       20.2
#> 4       1          MEP1     4                65          202       20.2
```
