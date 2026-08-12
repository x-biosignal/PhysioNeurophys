# Estimate the H-reflex recruitment threshold

Finds the lowest stimulation intensity whose mean H amplitude reaches an
absolute cutoff or a fraction of the observed Mmax.

## Usage

``` r
hReflexThreshold(recruitment, criterion = 0.05, scale = c("mmax", "absolute"))
```

## Arguments

- recruitment:

  An H-reflex recruitment table returned by
  [`hReflexRecruitment()`](https://x-biosignal.github.io/PhysioNeurophys/reference/hReflexRecruitment.md),
  or a data frame containing `intensity`, `h_p2p`, and `m_p2p`.

- criterion:

  Non-negative absolute amplitude, or a fraction in `[0, 1]` when
  `scale = "mmax"`.

- scale:

  Whether `criterion` is relative to Mmax or in native amplitude units.

## Value

A one-row data frame containing the first qualifying intensity, applied
cutoff, and scale.

## Examples

``` r
recruitment <- data.frame(
  intensity = c(10, 20, 30),
  h_p2p = c(0.1, 0.5, 1),
  m_p2p = c(0.2, 0.8, 2.5)
)
hReflexThreshold(recruitment, criterion = 0.2)
#>   threshold_intensity cutoff scale
#> 1                  20    0.5  mmax
```
