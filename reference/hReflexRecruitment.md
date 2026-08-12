# Measure an H-reflex recruitment series

Measures peak-to-peak direct motor (M) and Hoffmann-reflex (H) responses
in separate post-stimulus windows for each trial.

## Usage

``` r
hReflexRecruitment(
  x,
  intensity,
  m_window_ms = c(2, 12),
  h_window_ms = c(20, 45),
  epoch_start_ms = 0,
  channel = 1,
  assay_name = NULL
)
```

## Arguments

- x:

  A `PhysioExperiment` with time x channels x trials data.

- intensity:

  Finite stimulation intensity, one value per trial. Values remain in
  the caller's native stimulator unit.

- m_window_ms:

  Direct M-wave window relative to the stimulus, in milliseconds.

- h_window_ms:

  H-reflex window relative to the stimulus, in milliseconds.

- epoch_start_ms:

  Time of the first epoch sample relative to the stimulus.

- channel:

  One channel label or 1-based index.

- assay_name:

  Optional assay name.

## Value

A data frame with one row per trial and M/H peak-to-peak amplitudes and
fiducial samples.

## References

Palmieri RM, Ingersoll CD, Hoffman MA (2004). The Hoffmann reflex:
methodologic considerations and applications for use in sports medicine
and athletic training research. *Journal of Athletic Training*,
39:268-277.

## Examples

``` r
pe <- make_mep(n_trials = 3, intensities = c(30, 40, 50), seed = 1)
hReflexRecruitment(
  pe, c(10, 20, 30),
  m_window_ms = c(10, 18), h_window_ms = c(20, 60),
  epoch_start_ms = -20
)
#>   trial intensity      m_p2p m_pos_peak_sample m_neg_peak_sample     h_p2p
#> 1     1        10 0.07590877               171               155 0.1306161
#> 2     2        20 0.08112701               157               156 0.3022308
#> 3     3        30 0.06616223               165               151 1.0298831
#>   h_pos_peak_sample h_neg_peak_sample
#> 1               206               232
#> 2               216               242
#> 3               214               242
```
