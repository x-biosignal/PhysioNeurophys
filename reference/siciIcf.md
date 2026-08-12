# Summarize short-interval inhibition and intracortical facilitation

Summarize short-interval inhibition and intracortical facilitation

## Usage

``` r
siciIcf(paired, sici_range_ms = c(1, 5), icf_range_ms = c(7, 20))
```

## Arguments

- paired:

  A table returned by
  [`pairedPulse()`](https://x-biosignal.github.io/PhysioNeurophys/reference/pairedPulse.md)
  containing `isi_ms`, `n_conditioned`, and `ratio`.

- sici_range_ms:

  Closed SICI ISI range in milliseconds.

- icf_range_ms:

  Closed ICF ISI range in milliseconds.

## Value

A one-row data frame with conditioned-trial-weighted SICI and ICF
ratios, percentages, and trial counts.

## Examples

``` r
paired <- data.frame(
  isi_ms = c(2, 3, 10, 15),
  n_conditioned = rep(2L, 4),
  ratio = c(0.4, 0.5, 1.5, 1.4)
)
siciIcf(paired)
#>   sici_ratio sici_inhibition_pct n_sici icf_ratio icf_facilitation_pct n_icf
#> 1       0.45                  55      4      1.45                   45     4
```
