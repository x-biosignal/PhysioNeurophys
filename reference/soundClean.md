# Clean TMS-EEG sensor noise with SOUND

Implements simultaneous leave-one-sensor noise estimation followed by
the source-space Wiener reconstruction of the SOUND algorithm.

## Usage

``` r
soundClean(
  x,
  forward_model,
  lambda = 0.1,
  iterations = 10L,
  tolerance = 0.01,
  reference = "average",
  analysis_window_ms = NULL,
  epoch_start_ms = 0,
  assay_name = NULL,
  output_assay = "tms_sound",
  overwrite = FALSE
)
```

## Arguments

- x:

  A pulse-cleaned `PhysioExperiment`.

- forward_model:

  A list containing a finite channel-matched `leadfield`.

- lambda:

  Positive source-reconstruction regularization.

- iterations:

  Maximum number of complete simultaneous updates.

- tolerance:

  Non-negative maximum relative noise change for convergence.

- reference:

  Currently exactly `"average"`.

- analysis_window_ms:

  Optional closed noise-estimation interval. `NULL` uses all samples
  except recorded pulse/recharge gaps.

- epoch_start_ms:

  Time of the first epoch sample relative to the pulse.

- assay_name:

  Optional pulse-cleaned input assay.

- output_assay:

  Output assay name.

- overwrite:

  Whether to overwrite the selected input assay.

## Value

A modified `PhysioExperiment` with the fixed SOUND cleaning matrix,
complete noise history, and convergence diagnostics.

## Details

SOUND assumes diagonal sensor noise, an approximate channel-matched lead
field, and i.i.d. source amplitudes for its regularized reconstruction.
All channel noise estimates are updated simultaneously from the
preceding round, and one final cleaning matrix is applied to every
sample and trial. Pulse samples are excluded from noise estimation.
Regularization is not a significance threshold and is not optimized
against the reported TEP peaks.

## References

Mutanen TP, Metsomaa J, Liljander S, Ilmoniemi RJ (2018). Automatic and
robust noise suppression in EEG and MEG: The SOUND algorithm.
*NeuroImage*, 166:135-151.
[doi:10.1016/j.neuroimage.2017.10.021](https://doi.org/10.1016/j.neuroimage.2017.10.021)
