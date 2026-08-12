# Estimate H-reflex onset latency

Uses the baseline-standard-deviation threshold and persistent-run
detector from
[`mepLatency()`](https://x-biosignal.github.io/PhysioNeurophys/reference/mepLatency.md)
within an H-reflex-specific search window.

## Usage

``` r
hReflexLatency(
  x,
  h_window_ms = c(20, 45),
  epoch_start_ms = 0,
  baseline_window_ms = NULL,
  threshold_sd = 3,
  min_consecutive_ms = 1,
  channel = 1,
  assay_name = NULL
)
```

## Arguments

- x:

  A `PhysioExperiment` with time x channels x trials data.

- h_window_ms:

  H-reflex window relative to the stimulus, in milliseconds.

- epoch_start_ms:

  Time of the first epoch sample relative to the stimulus.

- baseline_window_ms:

  Optional baseline interval. The default uses all pre-stimulus samples.

- threshold_sd:

  Baseline standard-deviation multiplier.

- min_consecutive_ms:

  Required duration of a threshold crossing.

- channel:

  One channel label or 1-based index.

- assay_name:

  Optional assay name.

## Value

A data frame containing onset sample and stimulus-relative H-reflex
latency for each trial.

## Examples

``` r
pe <- make_mep(n_trials = 3, mep_latency_ms = 25, seed = 2)
hReflexLatency(pe, epoch_start_ms = -20)
#>   trial channel channel_label onset_sample h_latency_ms
#> 1     1       1          MEP1           NA           NA
#> 2     2       1          MEP1          227         25.2
#> 3     3       1          MEP1          227         25.2
```
