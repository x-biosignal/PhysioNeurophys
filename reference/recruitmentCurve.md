# Fit a Boltzmann MEP recruitment curve

Fits peak-to-peak MEP amplitude against stimulator intensity with a
three-parameter logistic curve.

## Usage

``` r
recruitmentCurve(
  x,
  channel = 1,
  response_window_ms = c(10, 60),
  epoch_start_ms = 0,
  intensity = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  An epoched MEP `PhysioExperiment`.

- channel:

  One channel label or 1-based index.

- response_window_ms:

  MEP response window in milliseconds.

- epoch_start_ms:

  Time of the first epoch sample relative to the pulse.

- intensity:

  Optional per-trial stimulation intensity in `%MSO`.

- assay_name:

  Optional assay name.

## Value

A `recruitment_curve` object with fitted parameters, predictions, fit
quality, and the underlying `nls` model.

## References

Devanne H, Lavoie BA, Capaday C (1997). Input-output properties and gain
changes in the human corticospinal pathway. *Experimental Brain
Research*, 114:329-338.
[doi:10.1007/PL00005641](https://doi.org/10.1007/PL00005641)

## Examples

``` r
pe <- make_mep(n_trials = 15, baseline_sd = 0.001, seed = 3)
recruitmentCurve(pe, epoch_start_ms = -20)
#> MEP recruitment curve
#>   MEP_max: 1.9974
#>   S50: 49.999 %MSO
#>   k: 5.0189 %MSO
#>   Maximum slope: 0.099495
#>   R-squared: 1
#>   Trials: 15; converged: TRUE
```
