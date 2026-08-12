# Clean TMS-EEG muscle artifact with SSP-SIR

Estimates one deterministic artifact subspace across trials and
reconstructs the retained activity through a channel-matched lead field.

## Usage

``` r
tmsSSPSIR(
  x,
  forward_model,
  artifact_window_ms = c(5, 50),
  n_artifact = NULL,
  variance_fraction = 0.9,
  lambda = 0.05,
  reference = "average",
  epoch_start_ms = 0,
  assay_name = NULL,
  output_assay = "tms_sspsir",
  overwrite = FALSE
)
```

## Arguments

- x:

  A pulse-cleaned `PhysioExperiment`.

- forward_model:

  A list containing a finite channel-matched `leadfield`.

- artifact_window_ms:

  Closed artifact-training interval.

- n_artifact:

  Optional exact number of artifact components.

- variance_fraction:

  Cumulative artifact variance target when `n_artifact` is `NULL`.

- lambda:

  Positive source-reconstruction regularization.

- reference:

  Currently exactly `"average"`.

- epoch_start_ms:

  Time of the first epoch sample relative to the pulse.

- assay_name:

  Optional pulse-cleaned input assay.

- output_assay:

  Output assay name.

- overwrite:

  Whether to overwrite the selected input assay.

## Value

A modified `PhysioExperiment` with fixed SSP-SIR reconstruction
metadata.

## Details

The pulse must already be removed or hardware absence must be explicitly
declared in `metadata(x)$tms_eeg$pulse_absent`. Both data and lead field
are average referenced with the same matrix. The returned data use the
fixed source-informed reconstruction `Lr G P`, not the plain SSP
projection `P`. SSP-SIR can attenuate neural signal when the artifact
and neural subspaces overlap; a lead field is an approximate physical
model, not proof that a retained deflection is cortical.

## References

Mutanen TP, Kukkonen M, Nieminen JO, et al. (2016). Recovering
TMS-evoked EEG responses masked by muscle artifacts. *NeuroImage*,
139:157-166.
[doi:10.1016/j.neuroimage.2016.05.028](https://doi.org/10.1016/j.neuroimage.2016.05.028)
