# Calculate the Hmax/Mmax ratio

Calculate the Hmax/Mmax ratio

## Usage

``` r
hMaxMMax(recruitment)
```

## Arguments

- recruitment:

  An H-reflex recruitment table returned by
  [`hReflexRecruitment()`](https://x-biosignal.github.io/PhysioNeurophys/reference/hReflexRecruitment.md),
  or a data frame containing `intensity`, `h_p2p`, and `m_p2p`.

## Value

A one-row data frame containing observed maxima, their lowest tied
intensities, and the dimensionless Hmax/Mmax ratio.

## Examples

``` r
recruitment <- data.frame(
  intensity = c(10, 20, 30),
  h_p2p = c(0.2, 1, 0.8),
  m_p2p = c(0.5, 1.5, 2)
)
hMaxMMax(recruitment)
#>   hmax hmax_intensity mmax mmax_intensity hmax_mmax_ratio
#> 1    1             20    2             30             0.5
```
