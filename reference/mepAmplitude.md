# Measure motor-evoked-potential amplitude

Calculates peak-to-peak MEP amplitude within a response window for every
selected channel and trial.

## Usage

``` r
mepAmplitude(
  x,
  response_window_ms = c(10, 60),
  epoch_start_ms = 0,
  channels = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A `PhysioExperiment` with time x channels x trials MEP data.

- response_window_ms:

  Two response-window limits relative to the TMS pulse, in milliseconds.

- epoch_start_ms:

  Time of the first epoch sample relative to the pulse.

- channels:

  Optional channel labels or 1-based indices.

- assay_name:

  Optional assay name.

## Value

A data frame with one row per selected channel and trial.

## References

Rossini PM, Burke D, Chen R, et al. (2015). Non-invasive electrical and
magnetic stimulation of the brain, spinal cord, roots and peripheral
nerves: basic principles and procedures for routine clinical and
research application. *Clinical Neurophysiology*, 126:1071-1107.
[doi:10.1016/j.clinph.2015.02.001](https://doi.org/10.1016/j.clinph.2015.02.001)

## Examples

``` r
pe <- make_mep(n_trials = 4, seed = 1)
mepAmplitude(pe, epoch_start_ms = -20)
#>   channel channel_label trial intensity_pct_mso       p2p pos_peak_sample
#> 1       1          MEP1     1          30.00000 0.1306161             206
#> 2       1          MEP1     2          43.33333 0.4681572             216
#> 3       1          MEP1     3          56.66667 1.5946432             214
#> 4       1          MEP1     4          70.00000 2.0055849             213
#>   neg_peak_sample pos_peak_time_ms neg_peak_time_ms
#> 1             232             21.0             26.2
#> 2             242             23.0             28.2
#> 3             238             22.6             27.4
#> 4             238             22.4             27.4
```
