# Store TMS pulse metadata

Stores one pulse row per trial without replacing other experiment
metadata. Intensity and resting motor threshold are expressed as percent
of maximum stimulator output (`%MSO`).

## Usage

``` r
setTMSpulses(
  x,
  intensity_pct_mso,
  target = NA_character_,
  coil = NA_character_,
  rmt_pct_mso = NA_real_,
  onset_sec = NA_real_
)
```

## Arguments

- x:

  A `PhysioExperiment` with a 2-D or epoched 3-D assay.

- intensity_pct_mso:

  Numeric stimulation intensity, one value per trial.

- target:

  Cortical target or target-muscle label.

- coil:

  Coil descriptor.

- rmt_pct_mso:

  Resting motor threshold in `%MSO`.

- onset_sec:

  Pulse onset in the source continuous recording, in seconds.

## Value

The modified `PhysioExperiment`.

## Examples

``` r
pe <- make_mep(n_trials = 3, intensities = c(40, 50, 60))
pe <- setTMSpulses(pe, c(40, 50, 60), target = "M1_FDI")
getTMSpulses(pe)
#>   pulse intensity_pct_mso target coil rmt_pct_mso onset_sec
#> 1     1                40 M1_FDI <NA>          NA        NA
#> 2     2                50 M1_FDI <NA>          NA        NA
#> 3     3                60 M1_FDI <NA>          NA        NA
```
