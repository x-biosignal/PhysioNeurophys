# Detect and summarize F-waves

Measures late F-wave peak-to-peak amplitude and onset latency after
peripheral stimulation, then reports persistence, latency dispersion,
and the mean F/M amplitude ratio.

## Usage

``` r
fWaveDetect(
  x,
  f_window_ms = c(20, 60),
  m_window_ms = c(2, 12),
  epoch_start_ms = 0,
  amplitude_threshold = 0.04,
  baseline_window_ms = NULL,
  threshold_sd = 3,
  min_consecutive_ms = 0.5,
  channel = 1,
  assay_name = NULL
)
```

## Arguments

- x:

  A `PhysioExperiment` with time x channels x trials data.

- f_window_ms:

  F-wave search window relative to the stimulus, in milliseconds.

- m_window_ms:

  Direct M-wave window relative to the stimulus, in milliseconds.

- epoch_start_ms:

  Time of the first epoch sample relative to the stimulus.

- amplitude_threshold:

  Minimum F-wave peak-to-peak amplitude in the assay's native unit.

- baseline_window_ms:

  Optional latency baseline interval. The default uses all pre-stimulus
  samples.

- threshold_sd:

  Baseline standard-deviation multiplier for onset.

- min_consecutive_ms:

  Required onset-threshold duration.

- channel:

  One channel label or 1-based index.

- assay_name:

  Optional assay name.

## Value

An `f_wave_result` list containing per-trial `waves`, a one-row
`summary`, and analysis `settings`.

## References

Fisher MA (2007). F-waves: physiology and clinical uses.
*TheScientificWorldJournal*, 7:144-160.
[doi:10.1100/tsw.2007.49](https://doi.org/10.1100/tsw.2007.49)

## Examples

``` r
pe <- make_mep(n_trials = 5, mep_latency_ms = 25, seed = 3)
fWaveDetect(
  pe, m_window_ms = c(10, 18), f_window_ms = c(20, 60),
  epoch_start_ms = -20
)
#> F-wave analysis
#>   Detected: 5/5; persistence: 100%
#>   Mean latency: 25.45 ms; chronodispersion: 1 ms
```
