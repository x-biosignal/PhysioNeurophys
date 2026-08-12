# Summarize paired-pulse MEP responses

Compares conditioned motor-evoked-potential amplitudes with the
aggregate unconditioned test response at each interstimulus interval
(ISI).

## Usage

``` r
pairedPulse(
  x,
  condition,
  isi_ms,
  response_window_ms = c(10, 60),
  epoch_start_ms = 0,
  channel = 1,
  assay_name = NULL,
  aggregate = c("mean", "median")
)
```

## Arguments

- x:

  An epoched `PhysioExperiment` containing MEP trials.

- condition:

  Per-trial labels, exactly `"test"` or `"conditioned"`.

- isi_ms:

  Per-trial interstimulus interval in milliseconds. Use `NA` for test
  trials and a positive finite value for conditioned trials.

- response_window_ms:

  MEP response window relative to the test pulse.

- epoch_start_ms:

  Time of the first epoch sample relative to the test pulse.

- channel:

  One channel label or 1-based index.

- assay_name:

  Optional assay name.

- aggregate:

  Whether to summarize amplitudes by the mean or median.

## Value

A data frame with one row per conditioned ISI and conditioned/test
amplitudes, their ratio, inhibition, and facilitation.

## References

Kujirai T, Caramia MD, Rothwell JC, et al. (1993). Corticocortical
inhibition in human motor cortex. *Journal of Physiology*, 471:501-519.
[doi:10.1113/jphysiol.1993.sp019912](https://doi.org/10.1113/jphysiol.1993.sp019912)

## Examples

``` r
pe <- make_mep(
  n_trials = 6,
  intensities = rep(50, 6),
  baseline_sd = 0,
  seed = 1
)
pairedPulse(
  pe,
  condition = c("test", "test", rep("conditioned", 4)),
  isi_ms = c(NA, NA, 2, 2, 10, 10),
  epoch_start_ms = -20
)
#>   isi_ms n_conditioned conditioned_amplitude n_test test_amplitude ratio
#> 1      2             2             0.9980267      2      0.9980267     1
#> 2     10             2             0.9980267      2      0.9980267     1
#>   inhibition_pct facilitation_pct
#> 1              0                0
#> 2              0                0
```
