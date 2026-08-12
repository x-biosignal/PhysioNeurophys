# Retrieve TMS pulse metadata

Retrieve TMS pulse metadata

## Usage

``` r
getTMSpulses(x)
```

## Arguments

- x:

  A `PhysioExperiment`.

## Value

A six-column data frame with one row per stored pulse, or a typed
zero-row table when TMS metadata are absent.

## Examples

``` r
getTMSpulses(make_mep(n_trials = 3))
#>   pulse intensity_pct_mso target coil rmt_pct_mso onset_sec
#> 1     1                30   <NA> <NA>          NA        NA
#> 2     2                50   <NA> <NA>          NA        NA
#> 3     3                70   <NA> <NA>          NA        NA
```
