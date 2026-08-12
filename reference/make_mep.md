# Simulate epoched TMS motor-evoked potentials

Generates biphasic MEP waveforms with known latency and a Boltzmann
amplitude-intensity relationship.

## Usage

``` r
make_mep(
  n_time = 600,
  n_channels = 1,
  n_trials = 12,
  sr = 5000,
  epoch_start_ms = -20,
  intensities = NULL,
  plateau = 2,
  s50 = 50,
  k = 5,
  mep_latency_ms = 20,
  mep_dur_ms = 10,
  baseline_sd = 0.02,
  seed = NULL
)
```

## Arguments

- n_time:

  Number of samples per epoch.

- n_channels:

  Number of MEP channels.

- n_trials:

  Number of TMS trials.

- sr:

  Sampling rate in Hz.

- epoch_start_ms:

  Time of the first sample relative to the TMS pulse.

- intensities:

  Per-trial stimulation intensity in `%MSO`.

- plateau, s50, k:

  Boltzmann recruitment parameters.

- mep_latency_ms:

  Planted MEP onset latency.

- mep_dur_ms:

  Duration of the one-cycle biphasic waveform.

- baseline_sd:

  Gaussian baseline-noise standard deviation.

- seed:

  Optional random seed.

## Value

An epoched `PhysioExperiment` with TMS pulse metadata.

## Examples

``` r
pe <- make_mep(n_trials = 5, seed = 1)
dim(SummarizedExperiment::assay(pe))
#> [1] 600   1   5
getTMSpulses(pe)
#>   pulse intensity_pct_mso target coil rmt_pct_mso onset_sec
#> 1     1                30   <NA> <NA>          NA        NA
#> 2     2                40   <NA> <NA>          NA        NA
#> 3     3                50   <NA> <NA>          NA        NA
#> 4     4                60   <NA> <NA>          NA        NA
#> 5     5                70   <NA> <NA>          NA        NA
```
