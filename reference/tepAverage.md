# Average TMS-evoked potentials and measure canonical peaks

Baseline-corrects every trial and channel before a mean or trimmed mean,
then finds signed peaks independently in each requested channel.

## Usage

``` r
tepAverage(
  x,
  components = c("N15", "P30", "N45", "P60", "N100", "P180"),
  windows_ms = NULL,
  channels = NULL,
  baseline_window_ms = c(-200, -20),
  epoch_start_ms = 0,
  average = c("mean", "trimmed"),
  trim = 0.2,
  assay_name = NULL
)
```

## Arguments

- x:

  A `PhysioExperiment` containing epoched TMS-EEG.

- components:

  Exact canonical component names in requested order.

- windows_ms:

  Optional named component x 2 override matrix.

- channels:

  Optional exact labels or 1-based indices.

- baseline_window_ms:

  Closed baseline interval.

- epoch_start_ms:

  Time of the first sample relative to the pulse.

- average:

  Exactly `"mean"` or `"trimmed"`.

- trim:

  Trim fraction in `[0, 0.5)`.

- assay_name:

  Optional assay name.

## Value

A `tms_tep` list with waveform, component, and settings tables.

## Details

Default signed peak windows are N15 (10–20 ms), P30 (20–40 ms), N45
(35–55 ms), P60 (50–75 ms), N100 (80–140 ms), and P180 (140–220 ms).
Exact ties use the earliest sample and a peak on a window endpoint is
marked as a boundary diagnostic. Canonical names define measurement
windows only. Auditory and somatosensory responses can remain, and no
component is evidence of a verified cortical source.

## References

Rogasch NC, Sullivan C, Thomson RH, et al. (2017). Analysing concurrent
transcranial magnetic stimulation and electroencephalographic data.
*NeuroImage*, 147:934-951.
[doi:10.1016/j.neuroimage.2016.10.031](https://doi.org/10.1016/j.neuroimage.2016.10.031)
