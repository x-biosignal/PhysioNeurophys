# Remove TMS pulse and recharge intervals

Replaces sample-exact closed pulse/recharge intervals by a linear or
cubic Hermite bridge, or blanks them to zero or a trial-specific
baseline. The selected input assay is preserved unless
`overwrite = TRUE`.

## Usage

``` r
removeTMSpulse(
  x,
  pulse_window_ms = c(-2, 10),
  recharge_windows_ms = NULL,
  method = c("cubic", "linear", "blank"),
  support_ms = 2,
  blank_value = c("baseline", "zero"),
  baseline_window_ms = c(-200, -20),
  epoch_start_ms = 0,
  assay_name = NULL,
  output_assay = "tms_pulse_clean",
  overwrite = FALSE
)
```

## Arguments

- x:

  A `PhysioExperiment` with time x channel or time x channel x trial
  data.

- pulse_window_ms:

  Increasing pulse-window limits in milliseconds.

- recharge_windows_ms:

  Optional increasing pair or n x 2 matrix.

- method:

  Exactly `"cubic"`, `"linear"`, or `"blank"`.

- support_ms:

  One positive derivative-support duration for cubic interpolation, in
  milliseconds.

- blank_value:

  Exactly `"baseline"` or `"zero"`.

- baseline_window_ms:

  Closed baseline interval in milliseconds.

- epoch_start_ms:

  Time of the first sample relative to the pulse.

- assay_name:

  Optional input assay name.

- output_assay:

  Output assay name.

- overwrite:

  Whether to overwrite the selected input assay.

## Value

A modified `PhysioExperiment` with an appended TMS-EEG audit step.

## Details

Closed intervals use
`round((time_ms - epoch_start_ms) / 1000 * sampling_rate) + 1`.
Overlapping or sample-adjacent pulse/recharge intervals are merged.
Linear and cubic methods require real samples on both sides and never
bridge across a second removed interval. Baseline settings do not affect
interpolation or zero blanking. Baseline blanking requires at least two
samples wholly outside every removed interval. The audit metadata
contains the complete altered mask, merged sample ranges, continuity
jumps, and robust baseline noise.

## References

Rogasch NC, Sullivan C, Thomson RH, et al. (2017). Analysing concurrent
transcranial magnetic stimulation and electroencephalographic data.
*NeuroImage*, 147:934-951.
[doi:10.1016/j.neuroimage.2016.10.031](https://doi.org/10.1016/j.neuroimage.2016.10.031)
